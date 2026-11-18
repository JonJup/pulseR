## Once-per-session warnings would fire in whichever test happens to run
## first and make WARN counts depend on file order. Mark them as shown;
## tests that check a warning reset its key themselves.
assign("local_row_order", TRUE, envir = .pulse_state)
