

% 定义参数
Lx = 100000; %m现实物理长度
Ly = 80000;  %m现实物理长度
Nx = 101; 
Ny = 101;
dx = Lx / (Nx - 1);  %m
dy = Ly / (Ny - 1);  %m
dt = 5;   %s
T =80000;    %s
Nt = ceil(T / dt);
g = 9.8; %m/s2
jingdu =100;
weidu =36.5;
jiaosudu=0.00007292;
f=2*jiaosudu*sin(weidu);%科氏力
k =1;%紊动运动粘性项系数
depth=20;  %平均水深


%{设置沉积参数
    rho_s = 2650;        % 颗粒密度 kg/m3 
    rho_f = 1000;        % 流体密度 kg/m3
    D = 0.0005;           % 颗粒直径 m
    nu = 1.004e-3;       % 流体动力粘度 N·s/㎡(牛顿秒/米方）  动力粘性系数
    settling = (rho_s - rho_f)*g*D^2/18/nu;    %斯托克方程，静水沉积速率
    ca=0.03;                    %糙率(曼宁系数)
    md=2500;                     %初始泥沙密度 kg/m3，一般砂岩
    %}目前使用的沉积流速
    vc=sqrt(2*g*D*(rho_s - rho_f)/rho_f);  %搬运临界流速
  % t=(1/ca)*(depth)^(1/6);%摩阻项系数（chezy系数）(已搬至下方进行精细计算)
   KX=0.2; %挟沙力系数
   MX=0.1; %挟沙力指数
   jiasu=1;  %地貌加速因子
   %}设置沉积参数

  %{设置风浪参数
 %风场1      %{设置风应力参数
        rho_air = 1.29; % 空气密度 kg/m3
       % cD=0.0005;     % 风曳力系数（已经搬至下方精细计算）
        windspeed=10;    % 平均风速
        jiaodu=-pi*1/4;      %风向与x轴正方向（正东向）的夹角
       %设置风应力参数}
  %风场2  
        windspeed2=5;        % 平均风速
        jiaodu2=pi*1/2;      %风向与x轴正方向（正东向）的夹角（南风）
  %风场3
        windspeed3=5;        % 平均风速
        jiaodu3=-pi*1/2;      %风向与x轴正方向（正东向）的夹角（北风）
 
  %风场4
        windspeed4=5;        % 平均风速
        jiaodu4=pi*3/4;      %风向与x轴正方向（正东向）的夹角（东南风）

       
       
       %{波浪辐射应力参数
     
     %计算波高  
         if 0<windspeed<=16.808
     waveH=-0.082+0.076*windspeed+0.011*windspeed^2;
        elseif 16.808<windspeed<40
     waveH=0.588+0.217*windspeed;
         end          
     waveE=(1/8)*rho_f*g*waveH^2;
     waveL=1;               %波长（瞎给的，需要斯托克斯波计算得出）
     waveK=1/waveL;          %波数
     waveAlpha=pi*3/2;      %波向角
     
     
     %波浪辐射应力参数}
 %}设置风浪参数
 
% 创建网格
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

% 初始化变量
laplacian_u = zeros(Nx, Ny);
laplacian_v = zeros(Nx, Ny);
damping_u = zeros(Nx, Ny);
damping_v = zeros(Nx, Ny);
velocity_modulus = zeros(Nx, Ny);
wind_u= zeros(Nx, Ny);
wind_v= zeros(Nx, Ny);
F=zeros(Nx, Ny);%吹程值
F2=zeros(Nx, Ny);
F3=zeros(Nx, Ny);
F4=zeros(Nx, Ny);
cD=zeros(Nx, Ny);
cD2=zeros(Nx, Ny);
cD3=zeros(Nx, Ny);
cD4=zeros(Nx, Ny);


waveN= zeros(Nx, Ny);
waveSxx =zeros(Nx, Ny);
waveSyy =zeros(Nx, Ny);
waveSxy =zeros(Nx, Ny);

s= zeros(Nx, Ny);       %单位网格泥沙浓度
c= zeros(Nx, Ny);       %单位网格沉积厚度
r= zeros(Nx, Ny);        %单位网格泥沙浓度改变量
e= zeros(Nx, Ny);        %单位网格更新后泥沙浓度
alls=zeros(Nx, Ny);     %单位网格不更新泥沙浓度
m=zeros(Nx, Ny);        %单位网格垂线挟沙力

chezy= zeros(Nx, Ny);%单位网格chezy系数


% 文件名（Topographic data）
filename ='C:\Users\ZENG\Desktop\daima\qinghaihu\qinghaihu.xyz';

% 打开文件
fileID = fopen(filename,'r');

% 检查文件是否打开成功
if fileID == -1
    error('Cannot open file: %s', filename);
end

% 读取数据
C = textscan(fileID, '%f %f %f');
fclose(fileID);

X1 = C{1};
Y1 = C{2};
Z1 = C{3};

% 创建新的网格
xq = linspace(min(X1), max(X1), Nx);
yq = linspace(min(Y1), max(Y1), Ny);
[Xq, Yq] = meshgrid(xq, yq);

% 插值
Zq = griddata(X1, Y1, Z1, Xq, Yq);

% 保存到矩阵d
d = Zq;
for i = 1:Nx
    for j = 1:Ny
        if d(i,j) < -30 || isnan(d(i,j))
            d(i,j) = 0;
        end
    end
end
% 使用中值滤波器进行噪音处理
filter_size = 3; % 值可以调整，根据具体情况来决定
d = medfilt2(d, [filter_size filter_size]);

d=d+1;
   for i = 1:Nx
    for j = 1:Ny
        if d(i,j) >=0.995
            d(i,j) = 1;
        end
    end
   end
% 初始化变量
z = ones(Nx, Ny);   %绘图水位
h = z - d;           
K= zeros(Nx, Ny);
u = zeros(Nx, Ny);
v = zeros(Nx, Ny);
b= zeros(Nx, Ny);   %流速标量
Q=zeros(Nx, Ny);   %流量标量
H=h+K;


%增加海心山{

d(52:54,49:51)=1;

%增加海心山}

%增加初始水深，以确保海心山内部计算错误不会影响周边{
for i = 1:Nx
    for j = 1:Ny
      h(52:54,49:51)=0;
    end
end


%增加初始水深，以确保海心山内部计算错误不会影响周边}




% 创建图像
figure;

% 设置流速场箭头长度和大小
arrow_length = 10;   
arrow_size = 1; 

Hn =H;
boolArray = true(Nx, Ny);
for i = 1:Nx
    for j = 1:Ny
        if h(i,j) > 0
            boolArray(i,j) = true;

        else
           boolArray(i,j) = false; 
        end
    end
end


%根据不同水深范围，给定床底糙率(相当于给出各网格精细的chezy系数){
hcaolv=h;
% 定义水深范围和对应的值<10m水深：糙率40；10-15m水深：糙率50；15-30m水深：糙率60
hrange= [10, 15, 30];
hvalues = [0.15, 0.1, 0.025];

%[0.1, 0.03, 0.025];
% 初始化一个与水深矩阵相同大小的结果矩阵
ca1 = zeros(size(hcaolv));

% 遍历不同水深范围
for i = 1:numel(hrange)
    hhrange = hrange(i);
    hvalue = hvalues(i);
    
    % 使用逻辑索引将符合条件的网格位置赋值
    ca1(hcaolv <= hhrange) = hvalue;
    
    % 重置水深矩阵中符合条件的网格，以便下一个范围
    hcaolv(hcaolv <= hhrange) = Inf; % 使用 Inf 表示已处理的网格
end


for i = 1:Nx
    for j = 1:Ny
        chezy(i,j)=(1/ca1(i,j))*(depth)^(1/6);%摩阻项系数（chezy系数）
    end
end
%根据不同水深范围，给定床底糙率}


%计算1风场吹程（
           
%给定吹程初始点（）
           sx = Nx-1;
           sy = 1;

           % 计算吹程

           for i = 1:Nx-1
              for j = 1:Ny-1
                   if boolArray(i, j)&& boolArray(i-1, j)&&boolArray(i+1, j)&&boolArray(i, j-1)&& boolArray(i, j+1)
                  F(i, j) = sqrt((i - sx)^2 * dy^2 + (j - sy)^2 * dx^2) ;
                  
                   end
              end
          end

%计算1风场吹程）

%计算2风场吹程（
           
%给定吹程初始点（）
          % sx2 = Nx-1;
          sx2 = 1; % 南风
           sy2 = Ny/2;

           % 计算吹程

           for i = 1:Nx-1
              for j = 1:Ny-1
                   if boolArray(i, j)&& boolArray(i-1, j)&&boolArray(i+1, j)&&boolArray(i, j-1)&& boolArray(i, j+1)
                  F2(i, j) = sqrt((i - sx2)^2 * dy^2 + (j - sy2)^2 * dx^2) ;
                   end
              end
          end

%计算2风场吹程）

%计算3风场吹程（
           
%给定吹程初始点（）
           sx3 = 1;
           sy3 = Ny/2;

           % 计算吹程

           for i = 1:Nx-1
              for j = 1:Ny-1
                   if boolArray(i, j)&& boolArray(i-1, j)&&boolArray(i+1, j)&&boolArray(i, j-1)&& boolArray(i, j+1)
                  F3(i, j) = sqrt((i - sx3)^2 * dy^2 + (j - sy3)^2 * dx^2) ;
                   end
              end
          end

%计算3风场吹程）
%计算4风场吹程（
           
%给定吹程初始点（）
           sx4 = 1;
           sy4 = Ny/2;

           % 计算吹程

           for i = 1:Nx-1
              for j = 1:Ny-1
                   if boolArray(i, j)&& boolArray(i-1, j)&&boolArray(i+1, j)&&boolArray(i, j-1)&& boolArray(i, j+1)
                  F4(i, j) = sqrt((i - sx4)^2 * dy^2 + (j - sy4)^2 * dx^2) ;
                   end
              end
          end

%计算4风场吹程）
     

boundaryPointsArray = identifyBoundaryPoints(boolArray);%确定边界网格
% 演化过程
for n = 1:Nt
   H=h+K;
   
    %布哈河入河口1         
   v(62:62,22:22) = -0.1;
   u(62:62,22:22) = 0.4;
    K(62:62,22:22) = 0 ;
    

    
    v(58:58,15:15) = -0.05;
   u(58:58,15:15) = -0.05;
    K(58:58,15:15) = 0 ;
    
    v(74:74,21:21) = 0.05;
   u(74:74,21:21) = 0.05;
    K(74:74,21:21) = 0 ;


    %沙柳河入河口
    v(91:91,56:56) = -0.05;
    u(91:91,56:56) = 0.1;
    K(91:91,56:56) = 0;
    
    
    h=1-d;    %重新更新水平面    
    H=h+K;    %重新更新水深
   for i = 1:Nx
    for j = 1:Ny
        if d(i,j) >=0.995
            d(i,j) = 1;
        end
    end
   end

    
    dn = d ;
    un = u;
    vn = v;
    bn=  b;
    Hn = H;
    Kn = K;
    zn = z;
    sn = s;
    cn = c;
    ctotal = c;%总高程改变量
    ctotaln = ctotal;
    rn = r;
    en = e;
    allsn=alls;
    mn=m;
    Qn=Q;
    
    laplacian_un = laplacian_u; 
    laplacian_vn = laplacian_v; 
    damping_un = damping_u;
    damping_vn = damping_v;
    velocity_modulusn = velocity_modulus;
    wind_un= wind_u;
    wind_vn= wind_v;


    Fn = F;%吹程值更新
    waveNn= waveN;
    waveSxxn =waveSxx;
    waveSyyn =waveSyy;
    waveSxyn =waveSxy;
    

    
    
   %流速场计算 
    for i = 2:Nx-1
        for j = 2:Ny-1
  

            
            if boolArray(i, j)&& boolArray(i-1, j)&&boolArray(i+1, j)&&boolArray(i, j-1)&& boolArray(i, j+1)
              
              %计算1风场风曳力系数
              if 5<=windspeed<=24.4
                cD(i,j)=(-0.65+0.384*log(((windspeed)^1.537)*(F(i,j)^0.177)*(Hn(i,j)^-0.026)))*10^-3;
                  
                 else if 0.5<=windspeed<5

                 cD(i,j)=(-0.65+0.384*log((5^1.537)*(F(i,j)^0.177)*(Hn(i,j)^-0.026)))/0.00086*(0.0283+(0.00513/windspeed))^2*10^-3;
                 end
                 end
         
             %计算2风场风曳力系数
              if 5<=windspeed2<=24.4
                cD2(i,j)=(-0.65+0.384*log(((windspeed2)^1.537)*(F2(i,j)^0.177)*(Hn(i,j)^-0.026)))*10^-3;
                  
                 else if 0.5<=windspeed2<5

                 cD2(i,j)=(-0.65+0.384*log((5^1.537)*(F2(i,j)^0.177)*(Hn(i,j)^-0.026)))/0.00086*(0.0283+(0.00513/windspeed2))^2*10^-3;
                 end
                 end
            %计算3风场风曳力系数
              if 5<=windspeed3<=24.4
                cD3(i,j)=(-0.65+0.384*log(((windspeed3)^1.537)*(F3(i,j)^0.177)*(Hn(i,j)^-0.026)))*10^-3;
                  
                 else if 0.5<=windspeed3<5

                 cD3(i,j)=(-0.65+0.384*log((5^1.537)*(F3(i,j)^0.177)*(Hn(i,j)^-0.026)))/0.00086*(0.0283+(0.00513/windspeed3))^2*10^-3;
                 end
                 end

            %计算4风场风曳力系数
              if 5<=windspeed4<=24.4
                cD4(i,j)=(-0.65+0.384*log(((windspeed4)^1.537)*(F4(i,j)^0.177)*(Hn(i,j)^-0.026)))*10^-3;
                  
                 else if 0.5<=windspeed4<5

                 cD4(i,j)=(-0.65+0.384*log((5^1.537)*(F4(i,j)^0.177)*(Hn(i,j)^-0.026)))/0.00086*(0.0283+(0.00513/windspeed4))^2*10^-3;
                 end
                 end


                
                laplacian_u(i,j) = (un(i,j+1) - 2*un(i,j) + un(i,j-1))/dx^2 + (un(i+1,j) - 2*un(i,j) + un(i-1,j))/dy^2;
                laplacian_v(i,j) = (vn(i,j+1) - 2*vn(i,j) + vn(i,j-1))/dx^2 + (vn(i+1,j) - 2*vn(i,j) + vn(i-1,j))/dy^2;

                velocity_modulus(i,j) = sqrt(un(i,j)^2 + vn(i,j)^2);
                damping_u(i,j) = g*(chezy(i,j)^-2)*velocity_modulus(i,j)*un(i,j)/Hn(i,j);
                damping_v(i,j) = g*(chezy(i,j)^-2)*velocity_modulus(i,j)*vn(i,j)/Hn(i,j);
                
                %（1风场）风应力项
                wind_u(i,j) = (rho_air*cD(i,j)*(windspeed)^2*cos(jiaodu))/rho_f/(Hn(i,j));
                wind_v(i,j) = (rho_air*cD(i,j)*(windspeed)^2*sin(jiaodu))/rho_f/(Hn(i,j));
                %（2风场）风应力项
                wind_u2(i,j) = (rho_air*cD2(i,j)*(windspeed2)^2*cos(jiaodu2))/rho_f/(Hn(i,j));
                wind_v2(i,j) = (rho_air*cD2(i,j)*(windspeed2)^2*sin(jiaodu2))/rho_f/(Hn(i,j));
               %（3风场）风应力项
                wind_u3(i,j) = (rho_air*cD3(i,j)*(windspeed3)^2*cos(jiaodu3))/rho_f/(Hn(i,j));
                wind_v3(i,j) = (rho_air*cD3(i,j)*(windspeed3)^2*sin(jiaodu3))/rho_f/(Hn(i,j));
               %（4风场）风应力项
                wind_u4(i,j) = (rho_air*cD4(i,j)*(windspeed4)^2*cos(jiaodu4))/rho_f/(Hn(i,j));
                wind_v4(i,j) = (rho_air*cD4(i,j)*(windspeed4)^2*sin(jiaodu4))/rho_f/(Hn(i,j));


               %控制风影响范围，从而山谷风、湖陆风不是全局风
               %控制西北风只能吹北部
               if F2(i,j)<40000
                   wind_u(i,j)=wind_u(i,j)/1;
                   wind_v(i,j)=wind_v(i,j)/1;
               end


                %波浪辐射应力
                waveN(i,j) =1/2*waveK*Hn(i,j)/sinh(2*waveK*Hn(i,j));
                waveSxx(i,j) =(waveE/2)*((2*waveN(i,j))*(1+(cos(waveAlpha))^2)-1);
                waveSyy(i,j) =(waveE/2)*((2*waveN(i,j))*(1+(sin(waveAlpha))^2)-1);
                waveSxy(i,j) =(waveE/2)*((waveN(i,j))*(sin(2*waveAlpha)));
                
                
                u(i,j) = un(i,j) - dt/dx*(0.5*(un(i,j+1)-un(i,j-1))*un(i,j)) ...
                                  - dt/dy*(0.5*(un(i+1,j)- un(i-1,j))*vn(i,j)) ...
                                  - g*dt/dx*(Kn(i,j+1) - Kn(i,j)) + f*dt*vn(i,j)...
                                  + dt*k*laplacian_u(i,j) - dt*damping_u(i,j)+dt*wind_u(i,j)+dt*wind_u2(i,j)+dt*wind_u3(i,j)+dt*wind_u4(i,j)...
                                   - dt/rho_f/Hn(i,j)*(   (waveSxx(i,j+1)-waveSxx(i,j-1))/(2*dx)+ (waveSxy(i+1,j)-waveSxy(i-1,j))/(2*dy)       ) ;

                v(i,j) = vn(i,j) - dt/dx*(0.5*(vn(i,j+1)-vn(i,j-1))*un(i,j)) ...
                                  - dt/dy*(0.5*(vn(i+1,j)-vn(i-1,j))*vn(i,j)) ...
                                  - g*dt/dy*(Kn(i+1,j) - Kn(i,j)) - f*dt*un(i,j)  ...
                                  + dt*k*laplacian_v(i,j) - dt*damping_v(i,j)+dt*wind_v(i,j)+dt*wind_v2(i,j)+dt*wind_v3(i,j)+dt*wind_v4(i,j)...
                           - dt/rho_f/Hn(i,j)*(   (waveSxy(i,j+1)-waveSxy(i,j-1))/(2*dx)+ (waveSyy(i+1,j)-waveSyy(i-1,j))/(2*dy)       ) ;
                              
                 %求解标量流速          
                 b(i,j)=sqrt(un(i,j)^2+vn(i,j)^2);
                  Q(i,j)=(dx*abs(un(i,j))+dy*abs(vn(i,j)));
            end
        end
    end
%  % 对湖泊边界流速进行处理
  %对湖泊边界流速进行处理{
       for i = 2:Nx-1
        for j = 2:Ny-1
           % if boolArray(i, j)&& boolArray(i-1, j)&&boolArray(i+1, j)&&boolArray(i, j-1)&& boolArray(i, j+1)
                laplacian_u(i,j) = (un(i,j+1) - 2*un(i,j) + un(i,j-1))/dx^2 + (un(i+1,j) - 2*un(i,j) + un(i-1,j))/dy^2;
                laplacian_v(i,j) = (vn(i,j+1) - 2*vn(i,j) + vn(i,j-1))/dx^2 + (vn(i+1,j) - 2*vn(i,j) + vn(i-1,j))/dy^2;

                velocity_modulus(i,j) = sqrt(un(i,j)^2 + vn(i,j)^2);
                damping_u(i,j) = g*(chezy(i,j)^-2)*velocity_modulus(i,j)*un(i,j)/Hn(i,j);
                damping_v(i,j) = g*(chezy(i,j)^-2)*velocity_modulus(i,j)*vn(i,j)/Hn(i,j);
                
                %风应力项
                wind_u(i,j) = (rho_air*cD(i,j)*(windspeed)^2*cos(jiaodu))/rho_f/(Hn(i,j));
                wind_v(i,j) = (rho_air*cD(i,j)*(windspeed)^2*sin(jiaodu))/rho_f/(Hn(i,j));
                
                %波浪辐射应力
                waveN(i,j) =1/2*waveK*Hn(i,j)/sinh(2*waveK*Hn(i,j));
                waveSxx(i,j) =(waveE/2)*((2*waveN(i,j))*(1+(cos(waveAlpha))^2)-1);
                waveSyy(i,j) =(waveE/2)*((2*waveN(i,j))*(1+(sin(waveAlpha))^2)-1);
                waveSxy(i,j) =(waveE/2)*((waveN(i,j))*(sin(2*waveAlpha)));
                
             if (boolArray(i, j)&&~boolArray(i, j-1))|(boolArray(i, j)&& ~boolArray(i, j+1))   %本身为湿，左或右方干网格
            %u的x方向作差的通量全为0，v的x方向通量为0
                 u(i,j) = 0 ;
                  

                v(i,j) = vn(i,j)   - dt/dy*(0.5*(vn(i+1,j)-vn(i-1,j))*vn(i,j)) ...
                                  - g*dt/dy*(Kn(i+1,j) - Kn(i,j)) - f*dt*un(i,j)  ...
                                  + dt*k*laplacian_v(i,j) - dt*damping_v(i,j)+dt*wind_v(i,j)...
                           - dt/rho_f/Hn(i,j)*(   (waveSxy(i,j+1)-waveSxy(i,j-1))/(2*dx)+ (waveSyy(i+1,j)-waveSyy(i-1,j))/(2*dy)       ) ;



            end
            
            if (boolArray(i, j)&&~boolArray(i+1, j))|(boolArray(i, j)&& ~boolArray(i-1, j))   %本身为湿，上或下方干网格
                %v的y方向作差的通量全为0，u的y方向通量为0  
               
                               
                 u(i,j) = un(i,j) - dt/dx*(0.5*(un(i,j+1)-un(i,j-1))*un(i,j)) ...
                                  - g*dt/dx*(Kn(i,j+1) - Kn(i,j)) + f*dt*vn(i,j)...
                                  + dt*k*laplacian_u(i,j) - dt*damping_u(i,j)+dt*wind_u(i,j)...
                                   - dt/rho_f/Hn(i,j)*(   (waveSxx(i,j+1)-waveSxx(i,j-1))/(2*dx)+ (waveSxy(i+1,j)-waveSxy(i-1,j))/(2*dy)       ) ;
                
                v(i,j) =0 ;             
            end
  
 
            if boolArray(i, j)&&((~boolArray(i+1, j)| ~boolArray(i-1, j))||(~boolArray(i, j-1)| ~boolArray(i, j+1))) 
            
                %上或下方有一个干网格且左或右也有一个
            u(i,j) =0;
           v(i,j) =0;
            
            end

           % end    
        end
    end

 

  %水位计算  
    for i = 2:Nx-1
        for j = 2:Ny-1
             if boolArray(i, j)&& boolArray(i-1, j)&&boolArray(i+1, j)&&boolArray(i, j-1)&& boolArray(i, j+1)
            K(i,j) = Kn(i,j) - dt/dx*(u(i,j)*Hn(i,j) - u(i,j-1)*Hn(i,j)) ...
                             - dt/dy*(v(i,j)*Hn(i,j) - v(i-1,j)*Hn(i,j));
            
            z(i,j) = K(i,j)+h(i,j) + d(i,j); 
           H(i,j)= K(i,j)+h(i,j); 
            
            if H(i,j)<0   %若出现水深为负值，则对其水位进行冻结赋值，给一个5cm的水位
               H(i,j)=0.05; 
           d(i,j)=-0.95;
            end  
             end
             
        end
    end

    %对湖泊边界水位进行处理{
       for i = 2:Nx-1
        for j = 2:Ny-1
           
            
            if (boolArray(i, j)&&~boolArray(i, j-1))|(boolArray(i, j)&& ~boolArray(i, j+1))   %本身为湿，左或右方干网格
             K(i,j) = Kn(i,j)- dt/dy*(v(i,j)*Hn(i,j) - v(i-1,j)*Hn(i,j));
            
           z(i,j) = K(i,j)+h(i,j) + d(i,j); 
           H(i,j)= K(i,j)+h(i,j); 
            end
            
               if (boolArray(i, j)&&~boolArray(i+1, j))|(boolArray(i, j)&& ~boolArray(i-1, j))   %本身为湿，上或下方干网格
             K(i,j) = Kn(i,j) - dt/dx*(u(i,j)*Hn(i,j) - u(i,j-1)*Hn(i,j));
            
           z(i,j) = K(i,j)+h(i,j) + d(i,j); 
           H(i,j)= K(i,j)+h(i,j); 
            end
  
              if boolArray(i, j)&&((~boolArray(i+1, j)| ~boolArray(i-1, j))||(~boolArray(i, j-1)| ~boolArray(i, j+1))) 
            
                %上或下方有一个干网格且左或右也有一个
             K(i,j) = Kn(i,j);
            
           z(i,j) = K(i,j)+h(i,j) + d(i,j); 
           H(i,j)= K(i,j)+h(i,j); 
 
            end
           if H(i,j)<0   %若出现水深为负值，则对其水位进行冻结赋值，给一个5cm的水位
               K(i,j)=0;
               H(i,j)=0.05; 
            end  
                
        end
    end
    
    
    %对湖泊边界水位进行处理}
    
     
% 根据 Hn 的幅度来选择颜色映射
colormap_values = (Hn > 0) .* abs(Hn) ; % 归一化 Hn 的幅度
colormap_jet = colormap(jet); % 使用 jet 颜色映射
colormap_jet(1, :) = [1, 1, 1]; % 将颜色映射的第一行设为白色
colormap(colormap_jet);
caxis([0 1]); % 设置颜色映射范围





    % 绘制水面高度
%  subplot(1, 3, 1);    
%   surf(X, Y, z, colormap_values); % 使用归一化的颜色映射值
 %   h1=colorbar;
%   set(get(h1,'Title'),'string','m');
 %   title(sprintf('Time Step: %d', n));
 %   xlabel('y');
 %   ylabel('x');
 %  zlabel('z');
 %   xlim([0 Lx]);
 %   ylim([0 Ly]);
 %  zlim([0 3]);
%   view(0, 60); 
  
   
   
    % 速度场（tsx）
   subplot(1, 3, 1);
    
    % 设置流线密度（density），这是一个控制流线密度的参数
   cla; % 清除当前坐标轴
    density =5; % 可根据需要进行调整
 % quiver(X, Y,  u, v, arrow_length, 'LineWidth', arrow_size);
 % quiver(X, Y, u, v, 'AutoScale', 'off', 'AutoScaleFactor', 2, 'LineWidth', 5);
   streamslice(X,Y,u,v,density);
  
    title(sprintf('矢量速度场 (Time Step: %d)', n));
    xlabel('y');
    ylabel('x');
   xlim([0 Lx]);
   ylim([0 Ly]);
    axis xy;

    
    
    % 矢量场（tsx）
   % subplot(1, 2, 2); 
   % quiver(X, Y,  u, v, arrow_length, 'LineWidth', arrow_size);
   % title(sprintf('矢量速度场 (Time Step: %d)', n));
  
  %  xlabel('x');
   % ylabel('y');
   % xlim([0 Lx]);
   % ylim([0 Ly]);
  
  
     
   %东西流速图
%  subplot(1, 3, 2);
%  pcolor(X,Y,u);
 %  shading interp;%插值绘制颜色
 %  h2=colorbar;
 %   set(get(h2,'Title'),'string','');
 % title(sprintf('东西速度场 (Time Step: %d)', n));
 %  xlabel('y');
 %   ylabel('x');
 %   xlim([0 Lx]);
 %  ylim([0 Ly]);
 %   axis xy;  

   %南北流速图
 % subplot(1, 3, 3);
 % pcolor(X,Y,v);
 %  shading interp;%插值绘制颜色
 %  h2=colorbar;
 %   set(get(h2,'Title'),'string','');
%  title(sprintf(['南北速度场 (Time Step: %d)'], n));
%   xlabel('y');
%    ylabel('x');
%    xlim([0 Lx]);
%   ylim([0 Ly]);
%    axis xy;  


   %速度云图
   subplot(1, 2, 2);
  pcolor(X,Y,b);
   shading interp;%插值绘制颜色
    h2=colorbar;
    set(get(h2,'Title'),'string','m/s');
   title(sprintf('速度场 (Time Step: %d)', n));
    xlabel('y');
   ylabel('x');
    xlim([0 Lx]);
    ylim([0 Ly]);
    axis xy;
  
 
   %曼宁系数图
 % subplot(2, 2, 3);
 % pcolor(X,Y,ca1);
%   shading interp;%插值绘制颜色
 %  h2=colorbar;
 %   set(get(h2,'Title'),'string','');
 % title(sprintf('速度场 (Time Step: %d)', n));
 %  xlabel('y');
 %   ylabel('x');
 %   xlim([0 Lx]);
 %  ylim([0 Ly]);
 %   axis xy;

   %风应力系数图
%  subplot(1, 2, 2);
 % pcolor(X,Y,cD2);
 %  shading interp;%插值绘制颜色
 %  h2=colorbar;
 %   set(get(h2,'Title'),'string','m');
%  title(sprintf('风应力10西北风 (Time Step: %d)', n));
%    xlabel('y');
%   ylabel('x');
%    xlim([0 Lx]);
%   ylim([0 Ly]);
%    axis xy;



   
 drawnow; 
 
 
end
function boundaryArray = identifyBoundaryPoints(existingBoolArray)
    [Nx, Ny] = size(existingBoolArray);
    boundaryArray = false(Nx, Ny); % 初始化为false
    
    for i = 2:Nx-1
        for j = 2:Ny-1
            neighbors = [existingBoolArray(i-1, j), existingBoolArray(i+1, j), existingBoolArray(i, j-1), existingBoolArray(i, j+1)];
            
            
             
            % 检查湿网格的情况
            if existingBoolArray(i, j) && ~all(neighbors) % 当前点是湿网格，但至少有一个邻居是干网格
                boundaryArray(i, j) = true;
            end
        end
    end
end
function isInteriorWet = isInteriorWetPoint(i, j, wetArray, boundaryArray)
    % 判断一个网格是否是内部湿网格
    isInteriorWet = wetArray(i, j) && ~boundaryArray(i, j);
end