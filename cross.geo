algebraic3d

solid box = orthobrick(-1, -1, -1; 1, 1, 1);
solid core1 = orthobrick(-0.05, -0.15, -0.05; 0.05, 0.15, 0.05 );
solid core2 = orthobrick(-0.15, -0.05, -0.05; 0.15, 0.05, 0.05 );

solid core = core1 or core2;
solid air = box and not core;

tlo core -col=[1,0,0];
tlo air -col=[0,0,1] -transparent;