# Stub for .consensus_partition(): the list shape cluster_consensus() reads.
mock_partition <- function(labels) {
        function(...) list(labels = labels, step = rep(1L, length(labels)),
                           salso = list(a = 1))
}

