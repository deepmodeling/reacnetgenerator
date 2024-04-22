#include <Python.h>

class C_Stack {
private:
  struct Node {
    int val;
    Node *prev;
  };
  Node *tail;

public:
  C_Stack();

  ~C_Stack();

  void push(int val);

  int pop();
};
