function  ls = getLs_cost1(data, ev1, ev2, ev12, ev21)
p = [];
[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Ov12',ev12);
[linerBot,linerUp,pointrBot,pointrUp,flagr] = getOuterTengents(data.Ov21',ev21);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
p1=cross(linelUp, linerUp);
p2=cross(linelBot,linerBot);
p1 = p1/p1(3);
p2 = p2/p2(3);
p=[p p1 p2];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Ov12',ev1);
[linerBot,linerUp,pointrBot,pointrUp,flagr] = getOuterTengents(data.Or',ev2);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
p1=cross(linelUp, linerUp);
p2=cross(linelBot,linerBot);
p1 = p1/p1(3);
p2 = p2/p2(3);
p=[p p1 p2];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Or',ev1);
[linerBot,linerUp,pointrBot,pointrUp,flagr] = getOuterTengents(data.Ov21',ev2);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
p1=cross(linelUp, linerUp);
p2=cross(linelBot,linerBot);
p1 = p1/p1(3);
p2 = p2/p2(3);
p=[p p1 p2];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Ov1',ev1);
[linerBot,linerUp,pointrBot,pointrUp,flagr] = getOuterTengents(data.Ov2',ev2);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
p1=cross(linelUp, linerUp);
p2=cross(linelBot,linerBot);
p1 = p1/p1(3);
p2 = p2/p2(3);
p=[p p1 p2];

num=size(p,2);
% v=LineRANSAC(p,num);
[u d v] = svd(p');
v=v(:,3);
ls= -sign(v(3))*v/norm(v(1:2));

fig=figure(5);
ShowPoly( data.Or, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
hold on
axis ij;
axis equal;
ShowPoly( data.Ov1, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15); 
ShowPoly( data.Ov2, 'FaceColor', MyPalette(3),  'EdgeColor',  MyPalette(3), 'FaceAlpha', 0.15);  
ShowPoly( data.Ov12, 'FaceColor', MyPalette(4),  'EdgeColor',  MyPalette(4), 'FaceAlpha', 0.15);  
ShowPoly( data.Ov21, 'FaceColor', MyPalette(5),  'EdgeColor',  MyPalette(5), 'FaceAlpha', 0.15);  
drawLine(ls, 'g-');
close(fig)






