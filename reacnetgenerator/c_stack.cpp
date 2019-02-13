extern "C"{
    #include "c_stack.h"
}

C_Stack::C_Stack() {
    tail = new Node;
    tail->prev = NULL;
    tail->val = -1;
};

C_Stack::~C_Stack() {
    Node *t;
    while(tail!=NULL){
        t=tail;
        tail=tail->prev;
        delete t;
    }
};

void C_Stack::push(int val) {
    Node* nt = new Node;
    nt->prev = tail;
    nt->val = val;
    tail = nt;
}

int C_Stack::pop() {
    Node* ot = tail;
    int val = tail->val;
    if (tail->prev != NULL) {
        tail = tail->prev;
        delete ot;
    }
    return val;
}