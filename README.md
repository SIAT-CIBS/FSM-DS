# FSM-DS
Dynamical Systems (DS) are a promising approach for robot motion modeling, which provide a flexible means to robot learning and control. Accuracy, stability and learning speed are the three major ingredients to learn robot motions from demonstrations with DS. Some approaches yield stable dynamical systems but potentially result in a poor reproduction performance, and some approaches yield good reproduction per- formance but are quite complex and time-consuming processes. Our focus in this paper is to address the accuracy-stability-speed issues simultaneously. We present a learning method named Fast and Stable Modeling for Dynamic Systems (FSM-DS), which is based on the Extreme Learning Machine to efficiently and accurately learn the parameters of the DS as will as to ensure the asymptotic stability at the target. The approach is validated through both 2-dimension tasks of learning handwriting motions and a set of robot experiments.

# [ ROS中常见坐标系定义及基本单位](http://blog.csdn.net/u010608582/article/details/52248115)  

为了方便开发者代码复用，ROS中统一定义了常见的坐标系（REP），所有的坐标系都是右手坐标系。
1. map
--固定的世界坐标系，z轴垂直向上。在map中表示的移动平台的pose是没有drift，没有累计误差的。而且，pos的变化是不连续的，
也就是说，pose是随时间间隔的跳动。
在典型应用中，系统中的定位模块接收到传感器观测值，不断的重复计算机器人的pose，这样就可以估计drift。但是在新的传感器信息
传递来时，会导致不连续的数据跳跃。
map参考系，作为long-term的全局参考是有用的，但是因为跳跃的存在，使得其对局部的行为的表达来说不是一个有用的参考。
2. odom
--固定的世界坐标系。在odom坐标系中，移动平台的pose是可以随时间无界的drift。这种drift使得odom在long-term的全局参考中不可用。
但是机器人的pose在odom中可以确保是连续的，意味着，变化是平滑的，无突变。
在典型应用中，odom的计算基于odometry信息源（包括了wheel odometry、visual odometry、IMU），在短时间内，基于odom的pose是
精确且连续的。
3. base_link
--固连在移动机器人上。base_link可以被配置在机器人的任何位置、任何角度。对不同的硬件而言，总是将其设置在对建模方便的位置角度上。
一般使用的坐标系设置如下：
x forward
y left
z up
相机坐标系：
z forward
x right
y down


4. 各参考系之间的关系。
map --> odom --> base_link
从odom到base_link的变换是基于任一odometry信息源，并且变换结果数据是由odometry发出。
从map到base_link的变换是由系统中的定位组件完成的。然而，变换结果不由定位组件发出。它仅仅接收odom->base_link的变换矩阵，然后发出
map->odom的变换矩阵
5. 数据单位
标准使用SI units
基本单位：
length: meter; mass: kilogram; time: second; current: ampere;
派生单位：
angle: radian; frequency: hertz; force: newton; power: watt; voltage: volt; temperature: celsius; magnetism: tesla;

注意RVIZ中使用的坐标系红是X轴，绿是Y轴，蓝是Z轴
```xml

<!-- 参数说明
static_transform_publisher x y z yaw pitch roll frame_id child_frame_id period_in_ms
 (yaw is rotation about Z, pitch is rotation about Y, and roll is rotation about X)
static_transform_publisher x y z qx qy qz qw frame_id child_frame_id  period_in_ms -->

<node pkg="tf" type="static_transform_publisher" name="link1_broadcaster" args="1 0 0 0 0 0 1 link1_parent link1 100" />
```
这里的yaw pitch roll是三个相对于父坐标系的旋转

参考：
http://www.ros.org/reps/rep-0103.html
http://www.ros.org/reps/rep-0105.html#id12
