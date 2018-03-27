function C = two_link_coriolis_matrix(dth1,dth2,l1,l2,m1,m2,th1,th2)
%TWO_LINK_CORIOLIS_MATRIX
%    C = TWO_LINK_CORIOLIS_MATRIX(DTH1,DTH2,L1,L2,M1,M2,TH1,TH2)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    10-Mar-2018 00:14:33

t2 = sin(th1);
t3 = cos(th1);
t4 = t3-1.0;
t5 = l1.*t3.*(1.0./2.0);
t14 = l1.*t4;
t6 = t5-t14;
t7 = sin(th2);
t8 = t3.*t7;
t9 = cos(th2);
t10 = t2.*t9;
t11 = t8+t10;
t12 = l1.*t3;
t13 = l1+l2;
t15 = l2.*(1.0./2.0);
t16 = l1+t15;
t17 = t2.*t7;
t21 = t3.*t9;
t22 = t17-t21;
t18 = t16.*t22;
t19 = t9-1.0;
t20 = t3.*t13.*t19;
t23 = -t17+t21;
t31 = t16.*t23;
t32 = t2.*t7.*t13;
t24 = t14+t20-t31-t32;
t25 = t23.*t24;
t26 = l1.*t11;
t27 = l1.*t2;
t28 = t3.*t7.*t13;
t29 = t2.*t13.*t19;
t48 = t11.*t16;
t30 = t27+t28+t29-t48;
t33 = t2.*t6;
t34 = l1.*t2.*t3.*(1.0./2.0);
t35 = t2.^2;
t36 = l1.*t35.*(1.0./2.0);
t37 = t3.*t6;
t38 = -t27+t33+t34;
t39 = t12+t20-t31-t32;
t40 = t11.*t24;
t41 = t9.*t16;
t60 = t13.*t19;
t42 = t41-t60;
t43 = t7.*t13;
t62 = t7.*t16;
t44 = t43-t62;
t45 = l1.*t23;
t50 = t23.*t39;
t46 = t25+t45-t50;
t65 = t11.*t39;
t47 = t26+t40-t65;
t49 = t11.*t30;
t53 = t23.*t30;
t51 = t26+t40-t53;
t52 = t25+t45+t49;
t54 = t2.*t9.*t13;
t55 = t28-t48+t54;
t71 = t3.*t9.*t13;
t56 = t31+t32-t71;
t57 = t3.^2;
t58 = l1.*t57.*(1.0./2.0);
t59 = -t12+t37+t58;
t61 = t9.*t13;
t63 = t7.*t44;
t70 = t9.*t42;
t64 = t61+t63-t70;
t66 = t7.*t42;
t67 = t41-t61;
t68 = t9.*t44;
t69 = -t43+t66+t68;
t72 = t11.*t56;
t73 = t23.*t55;
t74 = t26+t40-t53+t72+t73;
t75 = t23.*t56;
t80 = t11.*t55;
t76 = t25+t45+t49+t75-t80;
t77 = t9.*t67;
t78 = t61-t70+t77;
t79 = t7.*t67;
t81 = t43-t66+t79;
t82 = m2.*t64.*t74.*2.0;
t83 = m2.*t69.*t76.*2.0;
t84 = m2.*t52.*t81.*2.0;
t85 = m2.*t47.*t64.*2.0;
t86 = m2.*t46.*t69.*2.0;
t87 = m2.*t51.*t76.*2.0;
t88 = m2.*t64.*t74;
t89 = m2.*t69.*t76;
t90 = m2.*t52.*t81;
C = reshape([dth1.*(m1.*t38.*t59.*6.0+m2.*t46.*t51.*6.0-m2.*(t26-t11.*(t12+t18+t20-t2.*t7.*t13)+t11.*(t14+t18+t20-t2.*t7.*t13)).*(t25+t49-l1.*t22).*6.0+m1.*(t33+t34-l1.*t2).*(t12+t36-t3.*t6).*6.0).*(1.0./2.0)+dth2.*(t87+m1.*t38.*(t12+t36-t37).*2.0+m1.*t38.*t59.*2.0+m2.*t46.*t51.*2.0-m2.*t47.*t52.*2.0-m2.*t47.*t64-m2.*t46.*t69-m2.*t52.*t74.*2.0).*(1.0./2.0),dth1.*(t85+t86-t87+m2.*t52.*t74.*2.0).*(-1.0./2.0)-dth2.*(t82+t83+t84+m2.*t47.*t64+m2.*t46.*t69-m2.*t51.*t78.*2.0).*(1.0./2.0),dth1.*(t85+t86+t88+t89+t90-m2.*t51.*t78).*(-1.0./2.0)-dth2.*(t82+t83+t84-m2.*t51.*t78.*2.0).*(1.0./2.0),dth2.*(m2.*t64.*t81.*6.0+m2.*t69.*t78.*6.0).*(-1.0./2.0)-dth1.*(t88+t89+t90-m2.*t51.*t78+m2.*t64.*t81.*2.0+m2.*t69.*t78.*2.0).*(1.0./2.0)],[2,2]);
