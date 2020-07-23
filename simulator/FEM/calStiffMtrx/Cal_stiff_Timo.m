syms Tens Tors Bendy1 Bendy2 Bendy3 Bendy4 Bendz1 Bendz2 Bendz3 Bendz4
K = sym(zeros(12));
K(1,1) = Tens;
K(2,2) = Bendz1;
K(2,6) = Bendz2;
K(3,3) = Bendy1;
K(3,5) = -Bendy2;
K(4,4) = Tors;
K(5,5) = Bendy3;
K(6,6) = Bendz3;

K(1,7) = -Tens;
K(2,8) = -Bendz1;
K(2,12) = Bendz2;
K(3,9) = -Bendy1;
K(3,11) = -Bendy2;
K(4,10) = -Tors;
K(5,9) = Bendy2;
K(5,11) = Bendy4;
K(6,8) = -Bendz2;
K(6,12) = Bendz4;

K(7,7) = Tens;
K(8,8) = Bendz1;
K(8,12) = -Bendz2;
K(9,9) = Bendy1;
K(9,11) = Bendy2;
K(10,10) = Tors;
K(11,11) = Bendy3;
K(12,12) = Bendz3;

for i = 2:12
    for j = 1:(i-1)
        K(i,j) = K(j,i);
    end
end
Ke = K;

syms Q0 Q Qp0p L
Q0 = sym(zeros(1,12));
Q0(1) = Qp0p;
Q0(2:12) = sym('Qp%dp',[1,11]);
Q = sym(zeros(3));
QT = sym(zeros(3));
for i = 1:3
    for j = 1:3
        Q(i,j) = Q0(3*(j-1) + i);
        QT(j,i) = Q0(3*(j-1) + i);
    end
end
L = sym(zeros(12));
LT = sym(zeros(12));
for i = 1:4
    L((3*i-2):3*i,(3*i-2):3*i) = Q;
    LT((3*i-2):3*i,(3*i-2):3*i) = QT;
end
K = L*Ke*LT;

num = 0;
Matrix = sym(zeros(78,1));
for j = 1:12
    for i = j:-1:1
        num = num + 1;
        Matrix(num,1) = K(i,j);
    end
end
ccode(Matrix(:));