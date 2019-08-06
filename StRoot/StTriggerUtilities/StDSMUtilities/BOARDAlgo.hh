#ifndef BOARD_ALGO_HH
#define BOARD_ALGO_HH

struct Board;

struct BoardAlgo {
  virtual void operator()(Board& dsm) = 0;
  virtual ~BoardAlgo() {}
};

#endif	// BOARD_ALGO_HH

