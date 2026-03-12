#include <Rcpp.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <numeric>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

inline int to_type_id(const std::string& type) {
  if (type == "any") return 0;
  if (type == "start") return 1;
  if (type == "end") return 2;
  if (type == "within") return 3;
  if (type == "equal") return 4;
  Rcpp::stop("Unsupported `type`: %s", type);
}

inline bool strand_compatible(const int q_strand, const int s_strand, const bool ignore_strand) {
  if (ignore_strand) return true;
  if (q_strand == 0 || s_strand == 0) return true;
  return q_strand == s_strand;
}

inline bool interval_match(const int q_start,
                           const int q_end,
                           const int s_start,
                           const int s_end,
                           const int max_gap,
                           const int min_overlap,
                           const int type_id) {
  const int overlap_raw = std::min(q_end, s_end) - std::max(q_start, s_start) + 1;
  const int overlap_width = std::max(0, overlap_raw);

  int gap = -1;
  if (overlap_raw < 1) {
    if (s_start > q_end) {
      gap = s_start - q_end - 1;
    } else {
      gap = q_start - s_end - 1;
    }
  }

  if (max_gap < 0) {
    if (overlap_raw < 1) return false;
  } else {
    if (gap > max_gap) return false;
  }

  const int required_overlap = std::max(min_overlap, max_gap < 0 ? 1 : 0);
  if (overlap_width < required_overlap) return false;

  switch (type_id) {
  case 0:
    return true;
  case 1:
    if (max_gap < 0) return q_start == s_start;
    return std::abs(q_start - s_start) <= max_gap;
  case 2:
    if (max_gap < 0) return q_end == s_end;
    return std::abs(q_end - s_end) <= max_gap;
  case 3:
    if (!(q_start >= s_start && q_end <= s_end)) {
      return false;
    }
    if (max_gap <= 0) {
      return true;
    }
    return (s_end - s_start) - (q_end - q_start) <= max_gap;
  case 4:
    if (max_gap < 0) {
      return q_start == s_start && q_end == s_end;
    }
    return std::abs(q_start - s_start) <= max_gap &&
      std::abs(q_end - s_end) <= max_gap;
  default:
    return false;
  }
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_find_overlaps_indexed(
    const Rcpp::IntegerVector& q_start,
    const Rcpp::IntegerVector& q_end,
    const Rcpp::IntegerVector& q_seq,
    const Rcpp::IntegerVector& q_strand,
    const Rcpp::IntegerVector& s_start,
    const Rcpp::IntegerVector& s_end,
    const Rcpp::IntegerVector& s_seq,
    const Rcpp::IntegerVector& s_strand,
    const Rcpp::IntegerVector& s_original,
    const Rcpp::IntegerVector& partition_keys,
    const Rcpp::IntegerVector& partition_starts,
    const Rcpp::IntegerVector& partition_ends,
    const int max_gap,
    const int min_overlap,
    const std::string type,
    const bool ignore_strand,
    const int threads,
    const bool deterministic) {

  const int nq = q_start.size();
  const int ns = s_start.size();

  if (q_end.size() != nq || q_seq.size() != nq || q_strand.size() != nq) {
    Rcpp::stop("Query vectors must have equal length");
  }
  if (s_end.size() != ns || s_seq.size() != ns || s_strand.size() != ns || s_original.size() != ns) {
    Rcpp::stop("Subject vectors must have equal length");
  }
  if (partition_starts.size() != partition_ends.size() || partition_keys.size() != partition_starts.size()) {
    Rcpp::stop("Partition vectors must have equal length");
  }
  if (max_gap < -1) {
    Rcpp::stop("`max_gap` must be >= -1");
  }
  if (min_overlap < 0) {
    Rcpp::stop("`min_overlap` must be >= 0");
  }

  const int type_id = to_type_id(type);
  const int scan_gap = (max_gap < 0) ? 0 : max_gap;
  const int scan_pad = (max_gap < 0) ? 0 : 1;

  std::unordered_map<int, int> partition_by_key;
  partition_by_key.reserve(static_cast<std::size_t>(partition_keys.size() * 2));
  for (int i = 0; i < partition_keys.size(); ++i) {
    partition_by_key[partition_keys[i]] = i;
  }

  std::vector< std::vector<int> > query_partitions(static_cast<std::size_t>(partition_keys.size()));
  for (int qi = 0; qi < nq; ++qi) {
    const int seq_id = q_seq[qi];
    if (seq_id == 0) continue;
    const auto it = partition_by_key.find(seq_id);
    if (it == partition_by_key.end()) continue;
    query_partitions[static_cast<std::size_t>(it->second)].push_back(qi);
  }

  const std::size_t n_partitions = static_cast<std::size_t>(partition_keys.size());
  if (n_partitions == 0) {
    return Rcpp::List::create(
      Rcpp::Named("query_hits") = std::vector<int>(),
      Rcpp::Named("subject_hits") = std::vector<int>()
    );
  }

  struct Task {
    int part;
    int q_begin;
    int q_end;
  };

  const int thread_count = std::max(1, threads);
  const int min_chunk_size = 8192;
  const int target_chunks_per_thread = 4;

  std::vector< std::vector<int> > query_left_starts(n_partitions);
  std::vector<Task> tasks;
  tasks.reserve(n_partitions * 4);

  for (int part = 0; part < static_cast<int>(n_partitions); ++part) {
    std::vector<int>& q_idx = query_partitions[static_cast<std::size_t>(part)];
    if (q_idx.empty()) {
      continue;
    }

    std::sort(q_idx.begin(), q_idx.end(), [&](const int lhs, const int rhs) {
      if (q_start[lhs] != q_start[rhs]) return q_start[lhs] < q_start[rhs];
      if (q_end[lhs] != q_end[rhs]) return q_end[lhs] < q_end[rhs];
      return lhs < rhs;
    });

    const int s_begin = partition_starts[part] - 1;
    const int s_end_idx = partition_ends[part] - 1;
    std::vector<int>& left_starts = query_left_starts[static_cast<std::size_t>(part)];
    left_starts.resize(q_idx.size());

    if (s_begin < 0 || s_end_idx < s_begin) {
      std::fill(left_starts.begin(), left_starts.end(), s_end_idx + 1);
    } else {
      int left = s_begin;
      for (std::size_t q_pos = 0; q_pos < q_idx.size(); ++q_pos) {
        const int qi = q_idx[q_pos];
        const int ql = q_start[qi];
        while (left <= s_end_idx &&
               static_cast<long long>(s_end[left]) <
                 (static_cast<long long>(ql) - scan_gap - scan_pad)) {
          ++left;
        }
        left_starts[q_pos] = left;
      }
    }

    const int n_q = static_cast<int>(q_idx.size());
    const int desired_chunks = std::max(1, std::min(n_q, thread_count * target_chunks_per_thread));
    const int chunk_size = std::max(min_chunk_size, (n_q + desired_chunks - 1) / desired_chunks);
    for (int from = 0; from < n_q; from += chunk_size) {
      tasks.push_back(Task{
        part,
        from,
        std::min(n_q, from + chunk_size)
      });
    }
  }

  if (tasks.empty()) {
    return Rcpp::List::create(
      Rcpp::Named("query_hits") = std::vector<int>(),
      Rcpp::Named("subject_hits") = std::vector<int>()
    );
  }

  std::vector< std::vector<int> > out_query(tasks.size());
  std::vector< std::vector<int> > out_subject(tasks.size());

  auto process_task = [&](const std::size_t task_id) {
    const Task& task = tasks[task_id];
    const int part = task.part;
    const int s_begin = partition_starts[part] - 1;
    const int s_end_idx = partition_ends[part] - 1;

    if (s_begin < 0 || s_end_idx < s_begin) {
      return;
    }

    const std::vector<int>& q_idx = query_partitions[static_cast<std::size_t>(part)];
    const std::vector<int>& left_starts = query_left_starts[static_cast<std::size_t>(part)];
    std::vector<int>& q_hits = out_query[task_id];
    std::vector<int>& s_hits = out_subject[task_id];

    for (int q_pos = task.q_begin; q_pos < task.q_end; ++q_pos) {
      const int qi = q_idx[static_cast<std::size_t>(q_pos)];
      const int ql = q_start[qi];
      const int qr = q_end[qi];
      const int left = left_starts[static_cast<std::size_t>(q_pos)];
      if (left > s_end_idx) {
        continue;
      }

      for (int sj = left; sj <= s_end_idx; ++sj) {
        if (static_cast<long long>(s_start[sj]) >
              (static_cast<long long>(qr) + scan_gap + scan_pad)) {
          break;
        }

        if (!strand_compatible(q_strand[qi], s_strand[sj], ignore_strand)) {
          continue;
        }

        if (!interval_match(ql, qr, s_start[sj], s_end[sj], max_gap, min_overlap, type_id)) {
          continue;
        }

        q_hits.push_back(qi + 1);
        s_hits.push_back(s_original[sj]);
      }
    }
  };

  if (thread_count > 1 && tasks.size() > 1) {
    const int n_workers = std::min<int>(thread_count, static_cast<int>(tasks.size()));
    std::atomic<int> next_task(0);
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(n_workers));

    for (int worker = 0; worker < n_workers; ++worker) {
      workers.emplace_back([&]() {
        while (true) {
          const int task_id = next_task.fetch_add(1);
          if (task_id >= static_cast<int>(tasks.size())) {
            break;
          }
          process_task(static_cast<std::size_t>(task_id));
        }
      });
    }

    for (auto& worker : workers) {
      worker.join();
    }
  } else {
    for (std::size_t task_id = 0; task_id < tasks.size(); ++task_id) {
      process_task(task_id);
    }
  }

  std::size_t total_hits = 0;
  for (const auto& task_hits : out_query) {
    total_hits += task_hits.size();
  }

  std::vector<int> query_hits;
  std::vector<int> subject_hits;
  query_hits.reserve(total_hits);
  subject_hits.reserve(total_hits);

  for (std::size_t task_id = 0; task_id < tasks.size(); ++task_id) {
    query_hits.insert(query_hits.end(), out_query[task_id].begin(), out_query[task_id].end());
    subject_hits.insert(subject_hits.end(), out_subject[task_id].begin(), out_subject[task_id].end());
  }

  if (deterministic && query_hits.size() > 1) {
    std::vector<std::size_t> ord(query_hits.size());
    std::iota(ord.begin(), ord.end(), 0);

    std::sort(ord.begin(), ord.end(), [&](const std::size_t lhs, const std::size_t rhs) {
      if (query_hits[lhs] != query_hits[rhs]) return query_hits[lhs] < query_hits[rhs];
      if (subject_hits[lhs] != subject_hits[rhs]) return subject_hits[lhs] < subject_hits[rhs];
      return lhs < rhs;
    });

    std::vector<int> sorted_query(query_hits.size());
    std::vector<int> sorted_subject(subject_hits.size());
    for (std::size_t i = 0; i < ord.size(); ++i) {
      sorted_query[i] = query_hits[ord[i]];
      sorted_subject[i] = subject_hits[ord[i]];
    }

    query_hits.swap(sorted_query);
    subject_hits.swap(sorted_subject);
  }

  return Rcpp::List::create(
    Rcpp::Named("query_hits") = query_hits,
    Rcpp::Named("subject_hits") = subject_hits
  );
}

// [[Rcpp::export]]
Rcpp::IntegerVector cpp_count_overlaps_indexed(
    const Rcpp::IntegerVector& q_start,
    const Rcpp::IntegerVector& q_end,
    const Rcpp::IntegerVector& q_seq,
    const Rcpp::IntegerVector& q_strand,
    const Rcpp::IntegerVector& s_start,
    const Rcpp::IntegerVector& s_end,
    const Rcpp::IntegerVector& s_seq,
    const Rcpp::IntegerVector& s_strand,
    const Rcpp::IntegerVector& partition_keys,
    const Rcpp::IntegerVector& partition_starts,
    const Rcpp::IntegerVector& partition_ends,
    const int max_gap,
    const int min_overlap,
    const std::string type,
    const bool ignore_strand,
    const int threads) {

  const int nq = q_start.size();
  const int ns = s_start.size();

  if (q_end.size() != nq || q_seq.size() != nq || q_strand.size() != nq) {
    Rcpp::stop("Query vectors must have equal length");
  }
  if (s_end.size() != ns || s_seq.size() != ns || s_strand.size() != ns) {
    Rcpp::stop("Subject vectors must have equal length");
  }
  if (partition_starts.size() != partition_ends.size() || partition_keys.size() != partition_starts.size()) {
    Rcpp::stop("Partition vectors must have equal length");
  }
  if (max_gap < -1) {
    Rcpp::stop("`max_gap` must be >= -1");
  }
  if (min_overlap < 0) {
    Rcpp::stop("`min_overlap` must be >= 0");
  }

  const int type_id = to_type_id(type);
  const int scan_gap = (max_gap < 0) ? 0 : max_gap;
  const int scan_pad = (max_gap < 0) ? 0 : 1;
  const int thread_count = std::max(1, threads);
  const int min_chunk_size = 8192;
  const int target_chunks_per_thread = 4;

  std::unordered_map<int, int> partition_by_key;
  partition_by_key.reserve(static_cast<std::size_t>(partition_keys.size() * 2));
  for (int i = 0; i < partition_keys.size(); ++i) {
    partition_by_key[partition_keys[i]] = i;
  }

  std::vector< std::vector<int> > query_partitions(static_cast<std::size_t>(partition_keys.size()));
  for (int qi = 0; qi < nq; ++qi) {
    const int seq_id = q_seq[qi];
    if (seq_id == 0) continue;
    const auto it = partition_by_key.find(seq_id);
    if (it == partition_by_key.end()) continue;
    query_partitions[static_cast<std::size_t>(it->second)].push_back(qi);
  }

  const std::size_t n_partitions = static_cast<std::size_t>(partition_keys.size());
  std::vector<int> counts(static_cast<std::size_t>(nq), 0);
  if (n_partitions == 0 || nq == 0) {
    return Rcpp::wrap(counts);
  }

  struct Task {
    int part;
    int q_begin;
    int q_end;
  };

  std::vector< std::vector<int> > query_left_starts(n_partitions);
  std::vector<Task> tasks;
  tasks.reserve(n_partitions * 4);

  for (int part = 0; part < static_cast<int>(n_partitions); ++part) {
    std::vector<int>& q_idx = query_partitions[static_cast<std::size_t>(part)];
    if (q_idx.empty()) {
      continue;
    }

    std::sort(q_idx.begin(), q_idx.end(), [&](const int lhs, const int rhs) {
      if (q_start[lhs] != q_start[rhs]) return q_start[lhs] < q_start[rhs];
      if (q_end[lhs] != q_end[rhs]) return q_end[lhs] < q_end[rhs];
      return lhs < rhs;
    });

    const int s_begin = partition_starts[part] - 1;
    const int s_end_idx = partition_ends[part] - 1;
    std::vector<int>& left_starts = query_left_starts[static_cast<std::size_t>(part)];
    left_starts.resize(q_idx.size());

    if (s_begin < 0 || s_end_idx < s_begin) {
      std::fill(left_starts.begin(), left_starts.end(), s_end_idx + 1);
    } else {
      int left = s_begin;
      for (std::size_t q_pos = 0; q_pos < q_idx.size(); ++q_pos) {
        const int qi = q_idx[q_pos];
        const int ql = q_start[qi];
        while (left <= s_end_idx &&
               static_cast<long long>(s_end[left]) <
                 (static_cast<long long>(ql) - scan_gap - scan_pad)) {
          ++left;
        }
        left_starts[q_pos] = left;
      }
    }

    const int n_q = static_cast<int>(q_idx.size());
    const int desired_chunks = std::max(1, std::min(n_q, thread_count * target_chunks_per_thread));
    const int chunk_size = std::max(min_chunk_size, (n_q + desired_chunks - 1) / desired_chunks);
    for (int from = 0; from < n_q; from += chunk_size) {
      tasks.push_back(Task{
        part,
        from,
        std::min(n_q, from + chunk_size)
      });
    }
  }

  if (tasks.empty()) {
    return Rcpp::wrap(counts);
  }

  auto process_task = [&](const std::size_t task_id) {
    const Task& task = tasks[task_id];
    const int part = task.part;
    const int s_begin = partition_starts[part] - 1;
    const int s_end_idx = partition_ends[part] - 1;
    if (s_begin < 0 || s_end_idx < s_begin) {
      return;
    }

    const std::vector<int>& q_idx = query_partitions[static_cast<std::size_t>(part)];
    const std::vector<int>& left_starts = query_left_starts[static_cast<std::size_t>(part)];

    for (int q_pos = task.q_begin; q_pos < task.q_end; ++q_pos) {
      const int qi = q_idx[static_cast<std::size_t>(q_pos)];
      const int ql = q_start[qi];
      const int qr = q_end[qi];
      const int left = left_starts[static_cast<std::size_t>(q_pos)];
      if (left > s_end_idx) {
        counts[qi] = 0;
        continue;
      }

      int count = 0;
      for (int sj = left; sj <= s_end_idx; ++sj) {
        if (static_cast<long long>(s_start[sj]) >
              (static_cast<long long>(qr) + scan_gap + scan_pad)) {
          break;
        }

        if (!strand_compatible(q_strand[qi], s_strand[sj], ignore_strand)) {
          continue;
        }
        if (!interval_match(ql, qr, s_start[sj], s_end[sj], max_gap, min_overlap, type_id)) {
          continue;
        }

        ++count;
      }
      counts[qi] = count;
    }
  };

  if (thread_count > 1 && tasks.size() > 1) {
    const int n_workers = std::min<int>(thread_count, static_cast<int>(tasks.size()));
    std::atomic<int> next_task(0);
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(n_workers));

    for (int worker = 0; worker < n_workers; ++worker) {
      workers.emplace_back([&]() {
        while (true) {
          const int task_id = next_task.fetch_add(1);
          if (task_id >= static_cast<int>(tasks.size())) {
            break;
          }
          process_task(static_cast<std::size_t>(task_id));
        }
      });
    }

    for (auto& worker : workers) {
      worker.join();
    }
  } else {
    for (std::size_t task_id = 0; task_id < tasks.size(); ++task_id) {
      process_task(task_id);
    }
  }

  return Rcpp::wrap(counts);
}
