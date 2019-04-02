#include <string>
#include <list>
#include <vector>
#include <iostream>
#include <utility>
#include <math.h>
#include <ctime>

using namespace std;

#define PI 3.1415926
double angle(double angle) {
	return angle * PI / 180.0;
}


int virtual_location(int r_s, int e_n, pair<double, double> pos, double a, double b, double r_d) {
	vector<pair<double, double>> p_list; //��Ż���������λ����Ϣ
	p_list.push_back(pos);
	pair<double, double> v_pos;
	//cout << "weizhi1:" << pos.first << pos.second << endl;
	switch (e_n) {
	case 1:
	{
		if (abs(a - b) == 180) {
			for (int i = 1; i <= r_s / 2; i++) {
				v_pos.first = pos.first - i * r_d*cos(PI-angle(a));
				v_pos.second = pos.second - i * r_d*sin(PI-angle(a));
				cout <<"weizhi2:"<< v_pos.first << v_pos.second << endl;
				p_list.push_back(v_pos);
			}
			cout << "weizhi4:" << pos.first << pos.second << endl;
			for (int i = 1; i <= r_s / 2; i++) {
				v_pos.first = pos.first - i * r_d*cos(angle(b)-PI);
				v_pos.second = pos.second - i * r_d*sin(angle(b));
				cout << "weizhi3:" << v_pos.first <<"KK"<< v_pos.second << endl;
				p_list.push_back(v_pos);
			}
		}
		/*vector<pair<double, double>>::iterator iter;
		for (iter = p_list.begin(); iter != p_list.end(); iter++) {
			cout << iter->first << "/t" << iter->second << endl;
		}*/
		break;
	}
	case 2:
	{
		if (abs(a - b) != 180) {
			for (int i = 1; i <= r_s / 2; i++) {
				v_pos.first = pos.first - i * r_d*cos(PI-angle(a));
				v_pos.second = pos.second - i * r_d*sin(PI-angle(a));
				cout <<"weizhi4:"<< v_pos.first << v_pos.second << endl;
				p_list.push_back(v_pos);
			}

			for (int i = 1; i <= r_s / 2; i++) {
				v_pos.first = pos.first - i * r_d*sin(PI/2+angle(b));
				v_pos.second = pos.second - i * r_d*cos(PI/2+angle(b));
				cout << "weizhi5:" << v_pos.first <<"pp"<< v_pos.second << endl;
				p_list.push_back(v_pos);
			}
		}
		break;
	}
	case 3:
	{
		if (r_s % 3 != 0)
		{
			cout << "������������ζ��Σ�������β��ȶ���" << endl;
		}
		else{//����������Ӧ����3�ı�����������β��ȶ�����ʱ��r_s/e_nΪÿһ���߻����˵ĸ�����
			if (abs(a - b) == 60) {
				for (int i = 1; i <= r_s / e_n; i++) {
					cout << "weizhi4:" << pos.first << pos.second << endl;
					v_pos.first = pos.first - i * r_d*cos(PI-angle(a));
					v_pos.second = pos.second - i * r_d*sin(PI-angle(a));
					//cout << "weizhi4:" << v_pos.first << v_pos.second << endl;
					p_list.push_back(v_pos);
				}

				for (int i = 1; i <= r_s / e_n; i++) {
					v_pos.first = pos.first - i * r_d*cos(angle(b));
					v_pos.second = pos.second - i * r_d*sin(angle(b));
					//cout << "weizhi5:" << v_pos.first << v_pos.second << endl;
					p_list.push_back(v_pos);
				}

				vector<pair<double, double>>::iterator iter1;
				iter1 = p_list.end() - 1;
				pos.first = iter1->first;
				pos.second = iter1->second;
				//cout << "weizhi8:" << pos.first << pos.second << endl;

				for (int i = 1; i < r_s / e_n; i++) {
					v_pos.first = pos.first - i * r_d*sin(angle(b-a)-PI);
					v_pos.second = pos.second - i * r_d*cos(-angle(b-a)-PI);
					//cout << "weizhi9:" << v_pos.first << "hao" << v_pos.second << endl;
					p_list.push_back(v_pos);
				}
			}
		}
		break;
	}
	case 4:
	{
		if (r_s % 4 == 0) {//����������Ӧ����4�ı�����������β��ȶ�����ʱ��r_s/e_nΪÿһ���߻����˵ĸ�����
			for (int i = 1; i <= r_s / e_n; i++) {
				v_pos.first = pos.first - i * r_d*cos(PI-angle(a));
				v_pos.second = pos.second - i * r_d*sin(PI-angle(a));
				//cout <<"weizhi6:"<< v_pos.first << v_pos.second << endl;
				p_list.push_back(v_pos);
			}
			for (int i = 1; i <= r_s / e_n; i++) {
				v_pos.first = pos.first - i * r_d*sin(-angle(b));
				v_pos.second = pos.second - i * r_d*cos(-angle(b));
				//cout <<"weizhi7:"<< v_pos.first << v_pos.second << endl;
				p_list.push_back(v_pos);
			}

			vector<pair<double, double>>::iterator iter;
			iter = p_list.end()-1;
			pos.first = iter->first;
			pos.second = iter->second;
			//cout << "weizhi8:" << v_pos.first << v_pos.second << endl;
			
			for (int i = 1; i <= r_s / e_n; i++) {
				v_pos.first = pos.first - i * r_d*cos(PI-angle(a));
				v_pos.second = pos.second - i * r_d*sin(PI-angle(a));
				//cout << "weizhi9:" << v_pos.first << v_pos.second << endl;
				p_list.push_back(v_pos);
			}
			vector<pair<double, double>>::iterator iter1;
			iter1 = p_list.begin();
			for (int i = 0; i < r_s / e_n; i++) {
				iter1++;
			}
			pos.first = iter1->first;
			pos.second = iter1->second;
			//cout << "weizhi10:" << pos.first << pos.second << endl;

			for (int i = 1; i < r_s / e_n; i++) {
				v_pos.first = pos.first - i * r_d*sin(-angle(b));
				v_pos.second =pos.second - i * r_d*cos(-angle(b));
				//cout << "weizhi11:" << v_pos.first<<"hao" << v_pos.second << endl;
				p_list.push_back(v_pos);
			}
		}
		break;
	}
	default: 
	{
		if (r_s % e_n != 0) {
			cout << "�����γ��ȶ���" <<e_n<<"����"<< endl;
		}
		else
		{
			//ÿ�����ϻ����˵�����
			int num = r_s / e_n;
			v_pos.first = pos.first - num * r_d * cos(PI-angle(a));
			v_pos.second = pos.second -num * r_d * sin(PI-angle(a));
			cout << "jpy" << v_pos.first << "dd" << v_pos.second << endl;
			p_list.push_back(v_pos); //����������һ�����������λ�ã�����ȷ��Բ��

			//ȷ��Բ������
			double r;//�뾶
			r = r_d * r_s / (2 * e_n*sin(angle(180 / e_n)));
			cout <<"Բ�뾶Ϊ��" <<r << endl;

			//��֪Բ�����㡢�뾶��Բ������
			double x1, y1, x2, y2;
			double A, B, C, c1, c2;

			vector<pair<double, double>>::iterator iter;
			iter = p_list.begin();
			x1 = iter->first;
			y1 = iter->second;
			iter++;
			x2= iter->first;
			y2 = iter->second;
			cout << "x1=" << x1 << "y1=" << y1 << "x2=" << x2 << "y2=" << y2 << endl;

			//����Բ������
			double x_r, y_r;
			c1 = (x2*x2 - x1 * x1 + y2 * y2 - y1 * y1) /( 2 * (x2 - x1));
			c2 = (y2 - y1) / (x2 - x1);
			A = (c2 * c2 + 1);
			B = (2 * x1*c2 - 2 * c1*c2 - 2*y1);
			C = x1 * x1 - 2 * x1*c1 + c1 * c1 + y1* y1 - r * r;
			y_r = (-B + sqrt(B*B - 4 * A*C)) / (2 * A);
			x_r = c1 - c2 * y_r;
			cout << "Բ������Ϊ��x=" << x_r << "y=" << y_r << endl;

			//����������������Բ�������γɵļн�Ϊinc
			double inc = 360 / e_n;
			//��֪Բ�����꣬�캽�߻��������꣬�뾶�ͼнǣ���������̼��캽�����ꡣ




		}


	}
	}
	return 0;
}

int main() {
	int robot_size = 10;
	int edge_number = 5 ;
	pair<double, double> my_posoition = pair<double, double>(0, 0);
	double a = 30.0 ;
	double b = 90.0 ;
	double robot_distance = 5;
	virtual_location(robot_size, edge_number, my_posoition, a, b, robot_distance);
	return 0;
}
