function map = Set2(N)
% MatPlotLib 3.3 ��ɫ����
% ����:
% N   -  ����colormap���ȵ�������N>=0������Ϊ�գ���Ϊ��ǰͼ��colormap����
%
% ���:
% map -  Nx3��RGB��ɫ����
%
% Copyright  2020   Akun
% https://zhuanlan.zhihu.com/c_1074615528869531648

if nargin<1
	N = size(get(gcf,'colormap'),1);
else
	assert(isscalar(N)&&isreal(N),'First argument must be a real numeric scalar.')
	assert(fix(N)==N&&N>=0,'First argument must be a positive integer.')
end

C = [0.400000000000000,0.752941176470588,0.643137254901961;0.960784313725490,0.549019607843137,0.388235294117647;0.549019607843137,0.627450980392157,0.796078431372549;0.890196078431373,0.549019607843137,0.729411764705882;0.658823529411765,0.815686274509804,0.333333333333333;0.988235294117647,0.850980392156863,0.176470588235294;0.898039215686275,0.768627450980392,0.576470588235294;0.701960784313725,0.698039215686275,0.701960784313725];

map = C(1+mod(0:N-1,size(C,1)),:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����������            %%%
% ���ںţ������Ŀ����ճ� %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%