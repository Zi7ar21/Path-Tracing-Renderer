g++ -march=x86-64 -Og -s -std=c++11 -Wall -Wextra -Wshadow -fopenmp -I"./src/include" "./src/rng.cpp" "./src/main.cpp" -o "./OpenPT"