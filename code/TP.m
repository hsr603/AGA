function TP(D,U,Sc)
C=[D.Coord;D.Coord+Sc*U];e=D.Con(1,:);f=D.Con(2,:);
for i=1:6
    M=[C(i,e);C(i,f);repmat(NaN,size(e))];X(:,i)=M(:);    
end
plot3(X(:,1),X(:,2),X(:,3),'k',X(:,4),X(:,5),X(:,6),'m');axis('equal');
% plot3(X(:,1),X(:,2),X(:,3),'k');%����ǰ k��m������ɫ
% plot3(X(:,4),X(:,5),X(:,6),'m');%���κ�
axis('equal');

if D.Re(3,:)==1 %�����ƽ��ĸ˼� ���������Զ�ת���ӽ�
    view(2);
end

if 1 %*****************��ע�ڵ����
ccc=[1:size(D.Coord,2)];
text(D.Coord(1,:),D.Coord(2,:),D.Coord(3,:),num2str(ccc(:)));
clear ccc
%����Ҳ����
% for i=1:10
%    text(D.Coord(1,i),D.Coord(2,i),D.Coord(3,i),num2str(i))
% end
end


% grid on

%��������25�� �Ƿ�Ҳ�д�ϸ�Ĺ���
%�����Ϳ����ڻ�ͼ�����ֳ����˼��Ĵ�ϸ