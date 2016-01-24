#include "gravity.h"


void gravity(vector<double> &num,vector<double> &result) {
    const double pi = std::acos(-1.0);
    double angle = pi/num.size();
    Point temp(0.0, 0.0);
    vector<Point> point;
    vector<double> area;
    //double max_len =0.0;
    for (int i = 0; i < num.size(); ++i) {
        temp.x = std::cos(i*angle) * num[i];
        temp.y = std::sin(i*angle) * num[i];
        area.push_back(num[i]*num[(i+1)%num.size()]);
        point.push_back(temp);
     //   if(num[i]>max_len) max_len = num[i];
    }
    temp.x = -point[0].x;
    temp.y = 0;
    point.push_back(temp);
    double total_area = 0;
    Point grav(0.0, 0.0);
    temp = point[0];
    for (int i=0; i < area.size(); ++i) {
        point[i].x = (point[i].x + point[i+1].x)/3;
        point[i].y = (point[i].y + point[i+1].y)/3;
        total_area += area[i];
        grav.x += point[i].x*area[i];
        grav.y += point[i].y*area[i];
    }
    //for (int i=0; i<area.size();i++){
    //    result.push_back(area[i]/total_area);
    //}
    grav.x /= total_area;
    grav.y /= total_area;
    result.push_back(sqrt(grav.x*grav.x+grav.y*grav.y));
    result.push_back(std::atan2(grav.x,grav.y));
    result.push_back(total_area*sin(angle)/2);
}

