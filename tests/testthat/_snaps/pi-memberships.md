# Pi removes the small-region bias of normalised cohesion

    Code
      print(shown, row.names = FALSE)
    Output
       unit label mates_L mates_S cohesion_L cohesion_S old_U_L old_U_S Pi_L Pi_S
         40     L       5       5      0.128          1   0.114   0.886  0.5  0.5
         41     S       6       4      0.150          1   0.130   0.870  0.6  0.4

# memberships and summaries come from the repaired labels

    Code
      print(shown, row.names = FALSE)
    Output
       unit label_mocked label_repaired Pi1_if_unrepaired Pi1_returned certainty
          1            1              1             0.789        1.000     1.000
          2            1              1             0.770        0.996     0.996
          3            2              1             0.987        0.987     0.987
          4            1              1             0.707        0.979     0.979
          5            1              1             0.832        0.962     0.962
          6            1              1             0.872        0.958     0.958

# contiguity repair regrows freed units by Pi, not by cohesion

    Code
      print(shown, row.names = FALSE)
    Output
          score  to_L  to_S regrown_into
       cohesion 0.176 1.000            S
             Pi 0.545 0.455            L

