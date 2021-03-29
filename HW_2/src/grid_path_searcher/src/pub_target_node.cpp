#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/search/kdtree.h>
#include <pcl/search/impl/kdtree.hpp>

#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/PointCloud2.h>
#include <geometry_msgs/Vector3.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <Eigen/Eigen>
#include <math.h>
#include <random>

using namespace std;
using namespace Eigen;

ros::Publisher _target_pub;


int main(int argc, char** argv){
    ros::init(argc, argv, "pub_target_node");
    ros::NodeHandle nh("~");
    _target_pub = nh.advertise<geometry_msgs::PoseStamped>( "/goal", 1 );
    double x, y, z;
    nh.param("goal/x",        x,    4.64823);
    nh.param("goal/y",        y,    -1.80833);
    nh.param("goal/z",        z,    0.58);

    geometry_msgs::PoseStamped goal;
    int i = 0;
    ros::Rate rate(1);
    while (ros::ok())
    {
        goal.header.frame_id = std::string("world");
        goal.header.seq = 0;
        goal.header.stamp = ros::Time::now();
        goal.pose.position.x = x;
        goal.pose.position.y = y;
        goal.pose.position.z = z;
        _target_pub.publish(goal);
        i++;
        ros::spinOnce();
        rate.sleep();
    }
    return 0;
}