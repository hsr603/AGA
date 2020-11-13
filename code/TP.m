function TP(D,U,Sc)
C=[D.Coord;D.Coord+Sc*U];e=D.Con(1,:);f=D.Con(2,:);
for i=1:6
    M=[C(i,e);C(i,f);repmat(NaN,size(e))];X(:,i)=M(:);    
end
plot3(X(:,1),X(:,2),X(:,3),'k',X(:,4),X(:,5),X(:,6),'m');axis('equal');
% plot3(X(:,1),X(:,2),X(:,3),'k');%变形前 k和m控制颜色
% plot3(X(:,4),X(:,5),X(:,6),'m');%变形后
axis('equal');

if D.Re(3,:)==1 %如果是平面的杆件 ，这个则会自动转换视角
    view(2);
end

if 1 %*****************标注节点序号
ccc=[1:size(D.Coord,2)];
text(D.Coord(1,:),D.Coord(2,:),D.Coord(3,:),num2str(ccc(:)));
clear ccc
%这样也可以
% for i=1:10
%    text(D.Coord(1,i),D.Coord(2,i),D.Coord(3,i),num2str(i))
% end
end


% grid on

%后续看下25杆 是否也有粗细的规律
%这样就可以在绘图中体现出来杆件的粗细