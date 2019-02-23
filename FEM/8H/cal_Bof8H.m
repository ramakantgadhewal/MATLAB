function cal_Bof8H
%% % Construct GN
% % Construct symbolic xi, eta, zeta
% syms xi eta zeta GN;
% GN(1) = (1-eta)*(1-zeta)/8;
% GN(2) = (1+eta)*(1-zeta)/8;
% GN(3) = (1-eta)*(1+zeta)/8;
% GN(4) = (1+eta)*(1+zeta)/8;
% 
% GN(5) = (1+xi)*(1-zeta)/8;
% GN(6) = (1-xi)*(1-zeta)/8;
% GN(7) = (1+xi)*(1+zeta)/8;
% GN(8) = (1-xi)*(1+zeta)/8;
% 
% GN(9) = (1+xi)*(1-eta)/8;
% GN(10) = (1+xi)*(1+eta)/8;
% GN(11) = (1-xi)*(1+eta)/8;
% GN(12) = (1-xi)*(1-eta)/8;

%% % Construct GNe and Calculate Je
% % Construct symbolic GN and COORXYZ
% syms GN GNp0q;
% GN(1) = GNp0q;
% GN(2:12) = sym('GNp%dq',[1,11]);
% GNe = [GN(1)   GN(2)  -GN(2)  -GN(1)   GN(3)   GN(4)  -GN(4)  -GN(3);
%      -GN(5)   GN(5)   GN(6)  -GN(6)  -GN(7)   GN(7)   GN(8)  -GN(8);
%      -GN(9)  -GN(10)  -GN(11) -GN(12)  GN(9)   GN(10)   GN(11)  GN(12)];
% syms COORXYZp0q COORXYZ COORXYZe;
% COORXYZ(1) = COORXYZp0q;
% COORXYZ(2:24) = sym('COORXYZp%dq',[1 23]);
% COORXYZe = sym(zeros(8,3));
% for i = 1:8
%     for j = 1:3
%         NumFlag = 3*(i-1) + j;
%         COORXYZe(i,j) = COORXYZ(NumFlag);
%     end
% end
% % Calculate Je
% Je = GNe*COORXYZe;
% % Transform Je to a vector
% num = 0;
% J = sym(zeros(9,1));
% for j = 1:3
%     for i = 1:3
%         num = num + 1;
%         J(num,1) = Je(i,j);
%     end
% end
% J = ccode(J(:,1))

%% % Calculate invJe and detJ
% Construct symbolic J
syms J Je Jp0q;
J = sym(zeros(9,1));
J(1) = Jp0q;
J(2:9)= sym('Jp%dq',[1,8]);
Je = sym(zeros(3,3));
num = 0;
for j = 1:3
    for i = 1:3
        num = num + 1;
        Je(i,j) = J(num,1);
    end
end
% Calculate invJe and detJ
invJe = inv(Je);
detJ = det(Je)
% Transform invJe into a vector
num = 0;
invJ = sym(zeros(9,1));
for j = 1:3
    for i = 1:3
        num = num + 1;
        invJ(num,1) = invJe(i,j);
    end
end
invJ = ccode(invJ(:,1))

%% % Calculate kerBe
% % Construct symbolic invJ and GN
% syms invJ invJp0q GN GNp0q;
% invJ = sym(zeros(9,1));
% invJ(1) = invJp0q;
% invJ(2:9)= sym('invJp%dq',[1,8]);
% invJe = sym(zeros(3,3));
% num = 0;
% for j = 1:3
%     for i = 1:3
%         num = num + 1;
%         invJe(i,j) = invJ(num);
%     end
% end
% 
% GN = sym(zeros(12,1));
% GN(1) = GNp0q;
% GN(2:12) = sym('GNp%dq',[1,11]);
% GNe = [GN(1)   GN(2)  -GN(2)  -GN(1)   GN(3)   GN(4)  -GN(4)  -GN(3);
%       -GN(5)  GN(5)   GN(6)  -GN(6)  -GN(7)   GN(7)   GN(8)  -GN(8);
%       -GN(9)  -GN(10)  -GN(11) -GN(12)  GN(9)   GN(10)   GN(11)  GN(12)];
% 
% % Calculate kerBe
% kerBe = invJe*GNe;
% 
% % Transform kerBe into a vector
% kerB = sym(zeros(24,1));
% for j = 1:8
%     for i = 1:3
%         NumFlag = 3*(j-1) + i;
%         kerB(NumFlag,1) = kerBe(i,j);
%     end
% end
% kerB = ccode(kerB(:,1))

%% % Calculate Matrix
% % Construct symbolic kerB, D, detJ
% syms kerBp0q kerB D Dp0q detJ;
% kerB = sym(zeros(24,1));
% kerB(1) = kerBp0q;
% kerB(2:24)= sym('kerBp%dq',[1,23]);
% D = sym(zeros(3,1));
% D(1) = Dp0q;
% D(2:3)= sym('Dp%dq',[1,2]);
% % Calculate De, Be, BeT, Ke
% De = sym(zeros(6,6));
% De(1,1) = D(1);
% De(2,2) = D(1);
% De(3,3) = D(1);
% De(4,4) = D(3);
% De(5,5) = D(3);
% De(6,6) = D(3);
% De(1,2) = D(2);
% De(1,3) = D(2);
% De(2,1) = D(2);
% De(2,3) = D(2);
% De(3,1) = D(2);
% De(3,2) = D(2);
% Be= sym(zeros(6,24));
% BeT = sym(zeros(24,6));
% % node1
% Be(1:6,1:3) = [kerB(1)  0    0
%        0    kerB(2)  0
%        0    0   kerB(3)
%        kerB(2) kerB(1) 0
%        0 kerB(3) kerB(2)
%        kerB(3) 0 kerB(1)];
% % node2
% Be(1:6,4:6) = [kerB(4)  0    0
%        0    kerB(5)  0
%        0    0   kerB(6)
%        kerB(5) kerB(4) 0
%        0 kerB(6) kerB(5)
%        kerB(6) 0 kerB(4)];
% % node3
% Be(1:6,7:9) = [kerB(7)  0    0
%        0    kerB(8)  0
%        0    0   kerB(9)
%        kerB(8) kerB(7) 0
%        0 kerB(9) kerB(8)
%        kerB(9) 0 kerB(7)];
% % node4
% Be(1:6,10:12) = [kerB(10)  0    0
%        0    kerB(11)  0
%        0    0   kerB(12)
%        kerB(11) kerB(10) 0
%        0 kerB(12) kerB(11)
%        kerB(12) 0 kerB(10)];
% % node5
% Be(1:6,13:15) = [kerB(13)  0    0
%        0    kerB(14)  0
%        0    0   kerB(15)
%        kerB(14) kerB(13) 0
%        0 kerB(15) kerB(14)
%        kerB(15) 0 kerB(13)];
% % node6
% Be(1:6,16:18) = [kerB(16)  0    0
%        0    kerB(17)  0
%        0    0   kerB(18)
%        kerB(17) kerB(16) 0
%        0 kerB(18) kerB(17)
%        kerB(18) 0 kerB(16)];
% % node7
% Be(1:6,19:21) = [kerB(19)  0    0
%        0    kerB(20)  0
%        0    0   kerB(21)
%        kerB(20) kerB(19) 0
%        0 kerB(21) kerB(20)
%        kerB(21) 0 kerB(19)];
% % node8
% Be(1:6,22:24) = [kerB(22)  0    0
%        0    kerB(23)  0
%        0    0   kerB(24)
%        kerB(23) kerB(22) 0
%        0 kerB(24) kerB(23)
%        kerB(24) 0 kerB(22)];
% for i = 1:6
%     for j = 1:24
%         BeT(j,i) = Be(i,j);
%     end
% end
% Ke = BeT*De*Be*detJ;
% % Transform Ke to a vector called Matrix
% num = 0;
% Matrix = sym(zeros(300,1));
% for j = 1:24
%     for i = j:-1:1
%         num = num + 1;
%         Matrix(num,1) = Ke(i,j);
%     end
% end
% Matrix = ccode(Matrix(:))

%% % Calculate stressXYZ (i.e. stress on Gauss points)
% Construct symbolic kerB, D, Disp
syms kerB kerBAA0BB D DAA0BB Disp DispAA0BB
kerB = sym(zeros(24,1));
kerB(1) = kerBAA0BB;
kerB(2:24)= sym('kerBAA%dBB',[1,23]);
D = sym(zeros(3,1));
D(1) = DAA0BB;
D(2:3)= sym('DAA%dBB',[1,2]);
Disp = sym(zeros(24,1));
Disp(1) = DispAA0BB;
Disp(2:24)= sym('DispAA%dBB',[1,23]);
% Calculate De, Be, stressXYZe
De = sym(zeros(6,6));
De(1,1) = D(1);
De(2,2) = D(1);
De(3,3) = D(1);
De(4,4) = D(3);
De(5,5) = D(3);
De(6,6) = D(3);
De(1,2) = D(2);
De(1,3) = D(2);
De(2,1) = D(2);
De(2,3) = D(2);
De(3,1) = D(2);
De(3,2) = D(2);
Be= sym(zeros(6,24));
% node1
Be(1:6,1:3) = [kerB(1)  0    0
       0    kerB(2)  0
       0    0   kerB(3)
       kerB(2) kerB(1) 0
       0 kerB(3) kerB(2)
       kerB(3) 0 kerB(1)];
% node2
Be(1:6,4:6) = [kerB(4)  0    0
       0    kerB(5)  0
       0    0   kerB(6)
       kerB(5) kerB(4) 0
       0 kerB(6) kerB(5)
       kerB(6) 0 kerB(4)];
% node3
Be(1:6,7:9) = [kerB(7)  0    0
       0    kerB(8)  0
       0    0   kerB(9)
       kerB(8) kerB(7) 0
       0 kerB(9) kerB(8)
       kerB(9) 0 kerB(7)];
% node4
Be(1:6,10:12) = [kerB(10)  0    0
       0    kerB(11)  0
       0    0   kerB(12)
       kerB(11) kerB(10) 0
       0 kerB(12) kerB(11)
       kerB(12) 0 kerB(10)];
% node5
Be(1:6,13:15) = [kerB(13)  0    0
       0    kerB(14)  0
       0    0   kerB(15)
       kerB(14) kerB(13) 0
       0 kerB(15) kerB(14)
       kerB(15) 0 kerB(13)];
% node6
Be(1:6,16:18) = [kerB(16)  0    0
       0    kerB(17)  0
       0    0   kerB(18)
       kerB(17) kerB(16) 0
       0 kerB(18) kerB(17)
       kerB(18) 0 kerB(16)];
% node7
Be(1:6,19:21) = [kerB(19)  0    0
       0    kerB(20)  0
       0    0   kerB(21)
       kerB(20) kerB(19) 0
       0 kerB(21) kerB(20)
       kerB(21) 0 kerB(19)];
% node8
Be(1:6,22:24) = [kerB(22)  0    0
       0    kerB(23)  0
       0    0   kerB(24)
       kerB(23) kerB(22) 0
       0 kerB(24) kerB(23)
       kerB(24) 0 kerB(22)];
stressXYZe = De * Be * Disp;
% Transform stressXYZe to a vector stressXYZ
stressXYZe = ccode(stressXYZe)

%% % Calculate stressHex (stresses on the nodal points)
% % Construct symbolic interpo
% syms interpo interpoAA0BB stressXYZ stressXYZAA0BB;
% interpo = sym(zeros(4,1));
% interpo(1) = interpoAA0BB;
% interpo(2:4)= sym('interpoAA%dBB',[1,3]);
% stressXYZ = sym(zeros(8,1));
% stressXYZ(1) = stressXYZAA0BB;
% stressXYZ(2:8)= sym('stressXYZAA%dBB',[1,7]);
% % Calculate interpoe, recovery
% a = interpo(1);
% b = interpo(2);
% c = interpo(3);
% d = interpo(4);
% interpoe = [a, b, c, b, b, c, d, c,
%            b, a, b, c, c, b, c, d,
%            c, b, a, b, d, c, b, c,
%            b, c, b, a, c, d, c, b,
%            b, c, d, c, a, b, c, b,
%            c, b, c, d, b, a, b, c,
%            d, c, b, c, c, b, a, b,
%            c, d, c, b, b, c, b ,a];
% recoverye = interpoe * stressXYZ;
% % Transform recoverye to a vector recovery
% recovery = ccode(recoverye)
end