clc
clear all
close all
status=fclose('all');
% dbstop if error


if 1 %���ļ� д����
    d_Chrom=fopen('d_Chrom.txt','w+t');  %��¼�����Ʊ���
    d_Area=fopen('d_Area.txt','w+t');  %��¼Ӣ����� �й��ɵļ�¼
    d_ObjV=fopen('d_ObjV.txt','w+t');  %��¼Ŀ�꺯��
    d_F=fopen('d_F.txt','w+t');  %������Ӧ��
    d_U=fopen('d_U.txt','w+t');  %������λ��

    g_Chrom=fopen('G_Chrom.txt','w+t');
    g_Area=fopen('G_Area.txt','w+t');
    g_ObjV=fopen('G_ObjV.txt','w+t');
    g_F=fopen('G_F.txt','w+t');
    g_U=fopen('G_U.txt','w+t');
  
end

if 1  %��������
    global L p F_limit U_limit Area NVAR  my_rule_opst_num NIND
    
    D=Data10;

    Area=D.Area; p=D.p;U_limit=D.U_limit;F_limit=D.F_limit;
    T_rank=D.T; %�����ST�е�T�ظ�
    MAXGEN=99;
    GGAP=0.9;
    NIND=40;NVAR=size(T_rank,2);PRECI=1;
    GOOD_NUMSIZE=ceil(  NIND*(1-GGAP)); %��Ӣ������Ⱥ��С
    if GOOD_NUMSIZE<4 ;GOOD_NUMSIZE=4;  end %��Ϊranking������� �����ȥ�ĸ���Ҫ����3 ��˼����˵ ��Ϊ4
    BaseV=ones(1,NVAR)*max(size(Area));
    rule=ones(1,NVAR);
    my_rule_opst_num=rule;
    for i=1:size(D.Con,2) %����˼�����  ���ó�ȫ�ֱ��� �ߴ��Ż� �������Բ��� ��һ�ξͿ���
        H=D.Con(:,i);%�õ���i���˼��ı��
        C=D.Coord(:,H(2))-D.Coord(:,H(1)); %�˼������
        L(i)=norm(C);%norm(C)����������˼��ĳ���
    end
    Tpr_a=1-D.Re;f=find(Tpr_a); %��һ����Ϊ�˽�λ��ת�����õ��Ķ���
    clear C H Tpr_a
    NSA=zeros(MAXGEN+1,6);
end
tic
if 1 %��ʼ�� ����Chrom
    Initial_Chrom=crtbp( GOOD_NUMSIZE+NIND,BaseV);
    Initial_Chrom=Initial_Chrom+1;    %+1����Ϊ�����Ǹ����������ȱ��
    Initial_Chrom; %������A
    Initial_A=my10t25(Initial_Chrom,T_rank);
end

if 1 %��ʼ�� a)����Ŀ�꺯�� b)��Ӣ���ݿ� c)������ʼ��Ⱥ�ĸ���
    for i=1:size(Initial_A,1)  %����Ŀ�꺯��
        D.A=Initial_A(i,:)';
        [Tpr_F(i,:),d_u,Tpr_ObjV(i,:)]=my_objv(D);
        %        Tpr_U(i,:)=d_u(f)';
        Tpr_U(i,:)=reshape(d_u,1,3*size(d_u,2));
    end
    clear aa bb
    if 1  %��Ӣ��������
        Initial_FitnV=ranking(Tpr_ObjV(:,6),[2,1]) ;
        [~,num]=sort(Initial_FitnV,'descend');
        good=num(1:GOOD_NUMSIZE );
        G_Chrom (:,:,1)=Initial_Chrom(good,:);
        G_A(:,:,1)=Initial_A(good,:);
        G_ObjV(:,:,1)=Tpr_ObjV(good,:);
        G_F(:,:,1)=Tpr_F(good,:);
        G_U(:,:,1)=Tpr_U(good,:);
        GOOD=struct('Chrom',G_Chrom,'Area',G_A,'ObjV',G_ObjV,'F',G_F,'U',G_U);
    end
    if 1  %������Ϊ��Ӣ�������ǰ��ƶ��������ȺΪ������Ⱥ
        Initial_Chrom(good,:)=[];  Chrom(:,:,1)=Initial_Chrom;
        Initial_A(good,:)=[]; Unit_A(:,:,1)=Initial_A;
        Tpr_ObjV(good,:)=[];  ObjV(:,:,1)=Tpr_ObjV;
        Tpr_F(good,:)=[];     Unit_F(:,:,1)=Tpr_F;
        Tpr_U(good,:)=[];     Unit_U(:,:,1)=Tpr_U;
    end
    FitnV=ranking(ObjV(:,6),[2,1]) ;
    clear Initial_Chrom Initial_FitnV Initial_A
    
    U_trace(1,1)=min(ObjV(:,2,1));
    U_trace(1,2)=max(ObjV(:,2,1));
    U_trace(1,3)=mean(ObjV(:,2,1));
    U_trace(1,4)=std(ObjV(:,2,1));
    
end

if 1 %д������
    fprintf(d_Chrom,'��1��  \n '); fprintf(d_Area,'��1��  \n ');   fprintf(d_ObjV,'��1��  \n ');
    fprintf(d_F,'��1��  \n ');  fprintf(d_U,'��1��  \n ');
    for  i=1:NIND %DATA
        
        fprintf(d_Chrom,'#%d ',i) ;  fprintf(d_Chrom,'%d ',Chrom(i,:,1));   fprintf(d_Chrom,' \n ');
        fprintf(d_Area,'#%d ',i); fprintf(d_Area,'%.2f ',Area(Chrom(i,:,1))); fprintf(d_Area,' \n ');
        fprintf(d_ObjV,'#%d ',i);fprintf(d_ObjV,'%.2f ',ObjV(i,:,1)); fprintf(d_ObjV,' \n ');
        
        c=abs(Unit_F(i,:,1))/D.F_limit; index_F=find(c>1);
        d=abs(Unit_F(i,:,1))>D.F_limit;dd=double(d);dd(index_F)=c(index_F);
        fprintf(d_F,'#%d ',i);fprintf(d_F,'%d ',sum(d));fprintf(d_F,'(%d)',index_F);fprintf(d_F,' ');
        fprintf(d_F,'|%.2f| ',dd(index_F));  fprintf(d_F,'%.4f ',Unit_F(i,:,1));
        fprintf(d_F,' \n ');
        
        c=abs(Unit_U(i,:,1))/D.U_limit; index_U=find(c>1);
        d=abs(Unit_U(i,:,1))>D.U_limit;dd=double(d);dd(index_U)=c(index_U);
        fprintf(d_U,'#%d ',i);fprintf(d_U,'%d ',sum(d));fprintf(d_U,'(%d)',index_U);fprintf(d_U,' ');
        fprintf(d_U,'|%.2f| ',dd(index_U));  fprintf(d_U,'%.4f %.4f %.4f / ',Unit_U(i,:,1));
        fprintf(d_U,' \n ');
    end
    
    fprintf(g_Chrom,'��1��  \n '); fprintf(g_Area,'��1��  \n ');   fprintf(g_ObjV,'��1��  \n ');
    fprintf(g_F,'��1��  \n ');  fprintf(g_U,'��1��  \n ');
    for i=1:size(G_U,1) %GOOD
        fprintf(g_Chrom,'%d ',G_Chrom(i,:,1));   fprintf(g_Chrom,' \n ');
        fprintf(g_Area,'%.2f ',Area(G_Chrom(i,:,1))); fprintf(g_Area,' \n ');
        fprintf(g_ObjV,'%.2f ',G_ObjV(i,:,1)); fprintf(g_ObjV,' \n ');
        
        c=abs(G_F(i,:,1))/D.F_limit; index_F=find(c>1);
        d=abs(G_F(i,:,1))>D.F_limit;dd=double(d);dd(index_F)=c(index_F);
        fprintf(g_F,'%d ',sum(d));fprintf(g_F,'(%d)',index_F);fprintf(g_F,' ');
        fprintf(g_F,'|%.2f| ',dd(index_F));  fprintf(g_F,'%.4f ',G_F(i,:,1));
        fprintf(g_F,' \n ');
        
        c=abs(G_U(i,:,1))/D.U_limit; index_U=find(c>1);
        d=abs(G_U(i,:,1))>D.U_limit;dd=double(d);dd(index_U)=c(index_U);
        fprintf(g_U,'%d ',sum(d));fprintf(g_U,'(%d)',index_U);fprintf(g_U,' ');
        fprintf(g_U,'|%.2f| ',dd(index_U));  fprintf(g_U,'%.4f %.4f %.4f / ',G_U(i,:,1));
        fprintf(g_U,' \n ');
    end
    

end

for gen=1:MAXGEN
    %************ͨ��ѡ�� ���� ���� �õ��µĸ���
    Sel_Chrom=select('sus',Chrom(:,:,gen),FitnV,GGAP); %ѡ��
    Pc=0.8-0.6*gen/MAXGEN;  %����
    Sel_Chrom=xovmp(Sel_Chrom,Pc);Sel_Chrom=Sel_Chrom-1;
    %     Pm=0.2+0.4*gen/MAXGEN;  %����
    Pm=0.1+0.1*exp(-U_trace(gen,4)  )+0.2*gen/MAXGEN+0.1*exp(-std(  ObjV(:,4,gen)   )  );
    Sel_Chrom=mut(Sel_Chrom,Pm,BaseV);
    Son_Chrom=Sel_Chrom+1;  %�õ��Ӵ�����
    Son_A=my10t25(Son_Chrom, T_rank);
    
    for i=1:size(Son_A,1)   %�Ӵ�Ŀ�꺯��
        D.A=Son_A(i,:)';
        [~,~,Son_ObjV(i,:)]=my_objv(D);
    end
    
    if 1            %���ݵ�ѡ��ռ� �Ӵ��͸�����������
        Fat_ObjV=ObjV(:,:,gen);
        Chrom(:,:,gen+1)=reins(Chrom(:,:,gen),Son_Chrom,1,[1 0.1] ,Fat_ObjV(:,6),Son_ObjV(:,6) );
        Unit_A(:,:,gen+1)=my10t25( Chrom(:,:,gen+1), T_rank);
        for i=1:size(Unit_A,1)  %���¼���Ŀ�꺯��
            D.A=Unit_A(i,:,gen+1)';
            [Unit_F(i,:,gen+1),d_u,ObjV(i,:,gen+1)]=my_objv(D); %�Ӵ�Ŀ�꺯��
            %                Unit_U(i,:,gen+1)=d_u(f);
            Unit_U(i,:,gen+1)=reshape(d_u,1,3*size(d_u,2));
        end
    end


    U_trace(gen+1,1)=min(ObjV(:,2,gen+1));
    U_trace(gen+1,2)=max(ObjV(:,2,gen+1));
    U_trace(gen+1,3)=mean(ObjV(:,2,gen+1));
    U_trace(gen+1,4)=std(ObjV(:,2,gen+1));
    
    if 1 %���¾�Ӣ����
        FitnV=ranking(ObjV(:,6,gen+1),[2,1]) ;  %ĩβ��̭ �þ�Ӣ����������
        [~,num]=sort(FitnV,'ascend');
        if 1
            bad=num( +1:2*GOOD_NUMSIZE );
            Chrom(bad(1:GOOD_NUMSIZE),:,gen+1)=GOOD.Chrom(:,:,gen);
            Unit_A(bad(1:GOOD_NUMSIZE),:,gen+1)=GOOD.Area(:,:,gen);
            ObjV(bad(1:GOOD_NUMSIZE),:,gen+1)=GOOD.ObjV(:,:,gen);
            Unit_F(bad(1:GOOD_NUMSIZE),:,gen+1)=GOOD.F(:,:,gen);
            Unit_U(bad(1:GOOD_NUMSIZE),:,gen+1)=GOOD.U(:,:,gen);
        end
        
        FitnV=ranking(ObjV(:,6,gen+1),[2,1]) ;  %���¾�Ӣ����
        [~,num]=sort(FitnV,'descend');
        good=num(1: GOOD_NUMSIZE);
        G_Chrom (:,:,gen+1)=Chrom(good,:,gen+1);
        G_A(:,:,gen+1)=Unit_A(good,:,gen+1);
        G_ObjV(:,:,gen+1)=ObjV(good,:,gen+1);
        G_F(:,:,gen+1)=Unit_F(good,:,gen+1);G_U(:,:,gen+1)=Unit_U(good,:,gen+1);
        GOOD=struct('Chrom',G_Chrom,'Area',G_A,'ObjV',G_ObjV,'F',G_F,'U',G_U);
    end
    

    if 1 %д������
        fprintf(d_Chrom,'��%d��  \n ',gen+1); fprintf(d_Area,'��%d��  \n ',gen+1);fprintf(d_ObjV,'��%d��  \n ',gen+1);
        fprintf(d_F,'��%d��  \n ',gen+1);  fprintf(d_U,'��%d��  \n ',gen+1);
        for  i= 1:NIND %DATA
            fprintf(d_Chrom,'#%d ',i);  fprintf(d_Chrom,'%d ',Chrom(i,:,gen+1));   fprintf(d_Chrom,' \n ');
            fprintf(d_Area,'#%d ',i);fprintf(d_Area,'%.2f ',Area(Chrom(i,:,gen+1))); fprintf(d_Area,' \n ');
            fprintf(d_ObjV,'#%d ',i);fprintf(d_ObjV,'%.2f ',ObjV(i,:,gen+1)); fprintf(d_ObjV,' \n ');
            
            c=abs(Unit_F(i,:,gen+1))/D.F_limit; index_F=find(c> 1);
            d=abs(Unit_F(i,:,gen+1))>D.F_limit;dd=double(d);dd(index_F)=c(index_F);
            fprintf(d_F,'#%d ',i);fprintf(d_F,'%d ',sum(d));fprintf(d_F,'(%d)',index_F);fprintf(d_F,' ');
            fprintf(d_F,'|%.2f| ',dd(index_F));  fprintf(d_F,'%.4f ',Unit_F(i,:,gen+1));
            fprintf(d_F,' \n ');
            
            c=abs(Unit_U(i,:,gen+1))/D.U_limit; index_U=find(c> 1);
            d=abs(Unit_U(i,:,gen+1))>D.U_limit;dd=double(d);dd(index_U)=c(index_U);
            fprintf(d_U,'#%d ',i);fprintf(d_U,'%d ',sum(d));fprintf(d_U,'(%d)',index_U);fprintf(d_U,' ');
            fprintf(d_U,'|%.2f| ',dd(index_U));  fprintf(d_U,'%.4f %.4f %.4f / ',Unit_U(i,:,gen+1));
            fprintf(d_U,' \n ');
        end
        fprintf(g_Chrom,'��%d��  \n ',gen+1); fprintf(g_Area,'��%d��  \n ',gen+1);  fprintf(g_ObjV,'��%d��  \n ',gen+1);
        fprintf(g_F,'��%d��  \n ',gen+1);  fprintf(g_U,'��%d��  \n ',gen+1);
        for i= 1:size(G_U,1) %GOOD
            fprintf(g_Chrom,'%d ',G_Chrom(i,:,gen+1));   fprintf(g_Chrom,' \n ');
            fprintf(g_Area,'%.2f ',Area(G_Chrom(i,:,gen+1))); fprintf(g_Area,' \n ');
            fprintf(g_ObjV,'%.2f ',G_ObjV(i,:,gen+1)); fprintf(g_ObjV,' \n ');
            
            c=abs(G_F(i,:,gen+1))/D.F_limit; index_F=find(c> 1);
            d=abs(G_F(i,:,gen+1))>D.F_limit;dd=double(d);dd(index_F)=c(index_F);
            fprintf(g_F,'%d ',sum(d));fprintf(g_F,'(%d)',index_F);fprintf(g_F,' ');
            fprintf(g_F,'|%.2f| ',dd(index_F));  fprintf(g_F,'%.4f ',G_F(i,:,gen+1));
            fprintf(g_F,' \n ');
            
            c=abs(G_U(i,:,gen+1))/D.U_limit; index_U=find(c>1);
            d=abs(G_U(i,:,gen+1))>D.U_limit;dd=double(d);dd(index_U)=c(index_U);
            fprintf(g_U,'%d ',sum(d));fprintf(g_U,'(%d)',index_U);fprintf(g_U,' ');
            fprintf(g_U,'|%.2f| ',dd(index_U));  fprintf(g_U,'%.4f %.4f %.4f / ',G_U(i,:,gen+1));
            fprintf(g_U,' \n ');
        end
    end
    
end  %��ʼѭ���Ż�
toc
status=fclose('all');
disp('my_op.m end')




