clc
clear all
close all
status=fclose('all');
% dbstop if error


if 1 %打开文件 写数据
    d_Chrom=fopen('d_Chrom.txt','w+t');  %记录二进制编码
    d_Area=fopen('d_Area.txt','w+t');  %记录英寸面积 有规律的记录
    d_ObjV=fopen('d_ObjV.txt','w+t');  %记录目标函数
    d_F=fopen('d_F.txt','w+t');  %用来放应力
    d_U=fopen('d_U.txt','w+t');  %用来放位移

    g_Chrom=fopen('G_Chrom.txt','w+t');
    g_Area=fopen('G_Area.txt','w+t');
    g_ObjV=fopen('G_ObjV.txt','w+t');
    g_F=fopen('G_F.txt','w+t');
    g_U=fopen('G_U.txt','w+t');
  
end

if 1  %参数设置
    global L p F_limit U_limit Area NVAR  my_rule_opst_num NIND
    
    D=Data10;

    Area=D.Area; p=D.p;U_limit=D.U_limit;F_limit=D.F_limit;
    T_rank=D.T; %避免跟ST中的T重复
    MAXGEN=99;
    GGAP=0.9;
    NIND=40;NVAR=size(T_rank,2);PRECI=1;
    GOOD_NUMSIZE=ceil(  NIND*(1-GGAP)); %精英个体种群大小
    if GOOD_NUMSIZE<4 ;GOOD_NUMSIZE=4;  end %因为ranking这个函数 必须进去的个体要大于3 意思就是说 起步为4
    BaseV=ones(1,NVAR)*max(size(Area));
    rule=ones(1,NVAR);
    my_rule_opst_num=rule;
    for i=1:size(D.Con,2) %计算杆件长度  设置成全局变量 尺寸优化 拓扑属性不变 算一次就可以
        H=D.Con(:,i);%得到第i根杆件的编号
        C=D.Coord(:,H(2))-D.Coord(:,H(1)); %杆件坐标差
        L(i)=norm(C);%norm(C)是在求这根杆件的长度
    end
    Tpr_a=1-D.Re;f=find(Tpr_a); %这一步是为了将位移转换而用到的东西
    clear C H Tpr_a
    NSA=zeros(MAXGEN+1,6);
end
tic
if 1 %初始化 生成Chrom
    Initial_Chrom=crtbp( GOOD_NUMSIZE+NIND,BaseV);
    Initial_Chrom=Initial_Chrom+1;    %+1是因为上面那个函数本身的缺陷
    Initial_Chrom; %分组后的A
    Initial_A=my10t25(Initial_Chrom,T_rank);
end

if 1 %初始化 a)计算目标函数 b)精英数据库 c)修正初始种群的个数
    for i=1:size(Initial_A,1)  %计算目标函数
        D.A=Initial_A(i,:)';
        [Tpr_F(i,:),d_u,Tpr_ObjV(i,:)]=my_objv(D);
        %        Tpr_U(i,:)=d_u(f)';
        Tpr_U(i,:)=reshape(d_u,1,3*size(d_u,2));
    end
    clear aa bb
    if 1  %精英个体设置
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
    if 1  %修正因为精英个体而提前设计多出来的种群为正常种群
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

if 1 %写入数据
    fprintf(d_Chrom,'第1代  \n '); fprintf(d_Area,'第1代  \n ');   fprintf(d_ObjV,'第1代  \n ');
    fprintf(d_F,'第1代  \n ');  fprintf(d_U,'第1代  \n ');
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
    
    fprintf(g_Chrom,'第1代  \n '); fprintf(g_Area,'第1代  \n ');   fprintf(g_ObjV,'第1代  \n ');
    fprintf(g_F,'第1代  \n ');  fprintf(g_U,'第1代  \n ');
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
    %************通过选择 交叉 变异 得到新的个体
    Sel_Chrom=select('sus',Chrom(:,:,gen),FitnV,GGAP); %选择
    Pc=0.8-0.6*gen/MAXGEN;  %交叉
    Sel_Chrom=xovmp(Sel_Chrom,Pc);Sel_Chrom=Sel_Chrom-1;
    %     Pm=0.2+0.4*gen/MAXGEN;  %变异
    Pm=0.1+0.1*exp(-U_trace(gen,4)  )+0.2*gen/MAXGEN+0.1*exp(-std(  ObjV(:,4,gen)   )  );
    Sel_Chrom=mut(Sel_Chrom,Pm,BaseV);
    Son_Chrom=Sel_Chrom+1;  %得到子代编码
    Son_A=my10t25(Son_Chrom, T_rank);
    
    for i=1:size(Son_A,1)   %子代目标函数
        D.A=Son_A(i,:)';
        [~,~,Son_ObjV(i,:)]=my_objv(D);
    end
    
    if 1            %扩容的选择空间 子代和父代竞争生存
        Fat_ObjV=ObjV(:,:,gen);
        Chrom(:,:,gen+1)=reins(Chrom(:,:,gen),Son_Chrom,1,[1 0.1] ,Fat_ObjV(:,6),Son_ObjV(:,6) );
        Unit_A(:,:,gen+1)=my10t25( Chrom(:,:,gen+1), T_rank);
        for i=1:size(Unit_A,1)  %重新计算目标函数
            D.A=Unit_A(i,:,gen+1)';
            [Unit_F(i,:,gen+1),d_u,ObjV(i,:,gen+1)]=my_objv(D); %子代目标函数
            %                Unit_U(i,:,gen+1)=d_u(f);
            Unit_U(i,:,gen+1)=reshape(d_u,1,3*size(d_u,2));
        end
    end


    U_trace(gen+1,1)=min(ObjV(:,2,gen+1));
    U_trace(gen+1,2)=max(ObjV(:,2,gen+1));
    U_trace(gen+1,3)=mean(ObjV(:,2,gen+1));
    U_trace(gen+1,4)=std(ObjV(:,2,gen+1));
    
    if 1 %更新精英个体
        FitnV=ranking(ObjV(:,6,gen+1),[2,1]) ;  %末尾淘汰 用精英个体来代替
        [~,num]=sort(FitnV,'ascend');
        if 1
            bad=num( +1:2*GOOD_NUMSIZE );
            Chrom(bad(1:GOOD_NUMSIZE),:,gen+1)=GOOD.Chrom(:,:,gen);
            Unit_A(bad(1:GOOD_NUMSIZE),:,gen+1)=GOOD.Area(:,:,gen);
            ObjV(bad(1:GOOD_NUMSIZE),:,gen+1)=GOOD.ObjV(:,:,gen);
            Unit_F(bad(1:GOOD_NUMSIZE),:,gen+1)=GOOD.F(:,:,gen);
            Unit_U(bad(1:GOOD_NUMSIZE),:,gen+1)=GOOD.U(:,:,gen);
        end
        
        FitnV=ranking(ObjV(:,6,gen+1),[2,1]) ;  %更新精英个体
        [~,num]=sort(FitnV,'descend');
        good=num(1: GOOD_NUMSIZE);
        G_Chrom (:,:,gen+1)=Chrom(good,:,gen+1);
        G_A(:,:,gen+1)=Unit_A(good,:,gen+1);
        G_ObjV(:,:,gen+1)=ObjV(good,:,gen+1);
        G_F(:,:,gen+1)=Unit_F(good,:,gen+1);G_U(:,:,gen+1)=Unit_U(good,:,gen+1);
        GOOD=struct('Chrom',G_Chrom,'Area',G_A,'ObjV',G_ObjV,'F',G_F,'U',G_U);
    end
    

    if 1 %写入数据
        fprintf(d_Chrom,'第%d代  \n ',gen+1); fprintf(d_Area,'第%d代  \n ',gen+1);fprintf(d_ObjV,'第%d代  \n ',gen+1);
        fprintf(d_F,'第%d代  \n ',gen+1);  fprintf(d_U,'第%d代  \n ',gen+1);
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
        fprintf(g_Chrom,'第%d代  \n ',gen+1); fprintf(g_Area,'第%d代  \n ',gen+1);  fprintf(g_ObjV,'第%d代  \n ',gen+1);
        fprintf(g_F,'第%d代  \n ',gen+1);  fprintf(g_U,'第%d代  \n ',gen+1);
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
    
end  %开始循环优化
toc
status=fclose('all');
disp('my_op.m end')




