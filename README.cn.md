# Mono_AEC

#### 介绍
单声道音频回声消除(AEC)代码示例

使用自适应滤波器，根据两个信号输入（参考信号 x 和期望信号 d），找到合适的滤波器 w， 使得 d - x * w 最小。

其中，d = x * w0 + r，r为高斯白噪声，“ * ”表示卷积运算。

AEC函数代码改编自 Behrouz Farhang-Boroujeny 的书 Adaptive Filters:Theory and Applications，中文版请参考机械工业出版社：《自适应滤波器原理及Matlab仿真应用》

#### 软件架构
AECtest.m 作为演示主函数，调用5个AEC函数，分别是：
1.  VSNLMS:可变步长归一化最小均方算法
2.  VSNLMSNt:可变步长归一化最小均方牛顿法
3.  VSAPLMS:可变步长仿射投影算法
4.  VSNPFBLMS:可变步长归一化频域/快速分块最小均方算法
5.  SbLMS:子带最小均方算法

此外，PFBfilter及SbFilter分别在用于4.和5.中以训练好的参数对参考信号 x 进行滤波以逼近期望信号 d（由于这两种算法得到的自适应滤波器并非传统的卷积滤波器，不能简单使用filter函数进行滤波）
    
#### 安装教程

1.  将所有内容克隆到 Matlab 运行路径中
2.  运行 AECtest.m 


#### 使用说明

1.  运行 AECtest.m 并观察结果
2.  修改算法选项并观察结果
3.  修改输入音频/回声参数并观察结果


