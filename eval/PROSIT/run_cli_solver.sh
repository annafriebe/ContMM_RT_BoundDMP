# !/bin/bash

./cli_solver -R -T 200 -t 40 -q 6 -d 10 -m 6 -C ../../data/models/discr/em_distr6 -M ../../data/models/discr/tr_mat6.txt
./cli_solver -R -T 200 -t 50 -q 7 -d 8 -m 6 -C ../../data/models/discr/em_distr6 -M ../../data/models/discr/tr_mat6.txt
./cli_solver -R -T 200 -t 50 -q 8 -d 8 -m 6 -C ../../data/models/discr/em_distr6 -M ../../data/models/discr/tr_mat6.txt
