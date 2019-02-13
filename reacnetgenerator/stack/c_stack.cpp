extern "C"{
    #include "c_stack.h"
}

C_Stack::C_Stack() {
    tail = new Node;
    tail->prev = NULL;
    tail->val = NULL;
};

C_Stack::~C_Stack() {
    Node *t;
    while(tail!=NULL){
        t=tail;
        tail=tail->prev;
        delete t;
    }
};

PyObject* C_Stack::peek() {
    return tail->val;
}

void C_Stack::push(PyObject* val) {
    Node* nt = new Node;
    nt->prev = tail;
    nt->val = val;
    tail = nt;
}

PyObject* C_Stack::pop() {
    Node* ot = tail;
    PyObject* val = tail->val;
    if (tail->prev != NULL) {
        tail = tail->prev;
        delete ot;
    }
    return val;
}