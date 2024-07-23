#include "GuangXiSoftware.h"
#include "ui_GuangXiSoftware.h"

#include<Eigen/Dense>
#include <iostream>
#include <exception>
#include <cmath>
#include <limits>
#include <memory>

#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

#include <stdlib.h>
#include <string>
#include <QFont>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QProgressBar>
#include <QPushButton>
#include <QGridLayout>
#include <QProgressDialog>
#include <qlogging.h>
#include <QFont>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QProgressBar>
#include <QPushButton>
#include <QGridLayout>
#include <QProgressDialog>

#include "las_io.h"
#include "point_types.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <QChart>
#include <QScatterseries>
#include <QPolarchart>

#define e 2.718281828459
#define PI 3.1415926535
using namespace std;
using namespace Eigen;
using PointT = wl::LASPoint;
constexpr double DEG_TO_RAD_LOCAL = 3.1415926535897932 / 180.0;

struct PointPOS
{
	double time, latitude, longitude, altitude, x_velocity, y_velocity, z_velocity, roll, pitch, platform_heading;
	double  wander_angle, x_a, y_a, z_a, x_ar, y_ar, z_ar;
};
vector<PointPOS> points;

struct RXPPOS
{
	double X, Y, Z;
	double range, theta;
	double heading, pitch, roll;
	double wander;
	double latitude, longitude, altitude;

	double UTM_E;   //投影到 UTM 坐标系的横坐标
	double UTM_N;   //投影到 UTM 坐标系的纵坐标
	double UTM_H;   // 点云的高度
}q;
vector<RXPPOS> pointsRXPPOS;
double time_index = 0;

// 椭球参数
double A = 6378137.0;
double B = 6356752.3141;
double e2 = 0.0066943799013;
double N;
double GPSX, GPSY, GPSZ;

// 一个文件中 5500 帧（帧头有效），每帧 8 个通道回波
double rad2deg(double rad)
{
	rad = rad * 180 / PI;
	return rad;
}
double deg2rad(double deg)
{
	deg = deg * PI / 180;
	return deg;
}

int framenum = 0;     //5500   34
// 每帧长度
int framelength = 0;  //32768  480492
// 共计 30k 个点组成回波
int wavepoints = 0;   //2000   30000
// 主波回波起始字节位置
int mainwaveIndex = 0;//144    138
int waveIndex = 0;    //8204   60183

// 找数组最大值及其索引
void max(double* arr, double* max, int* maxindex)
{
	double tempMax = arr[0];
	int tempMaxIndex = 0;
	for (int i = 0; i < sizeof(arr) / sizeof(double); i++)
	{
		if (arr[i] > tempMax)
		{
			tempMax = arr[i];
			tempMaxIndex = i;
		}
	}
	*max = tempMax;
	*maxindex = tempMaxIndex;
}

// WGS-84 转经纬高 BLH 坐标系
void xyz2BLH(double x, double y, double z, double& longi, double& latit, double& Height)
{
	double p = sqrt(x * x + y * y);
	double angle = atan(z * A / (p * B));
    longi = atan(y / x);

/*	// 迭代公式
	double latit0 = atan(z / ((1 - e2) * p));
	double NN1 = A / sqrt(1 - e2 * sin(latit0) * sin(latit0));
	Height = p / cos(latit0) - NN1;
	latit = atan(z / ((1 - e2 * NN1 / (NN1 + Height) * p)));
*/	
    // 近似公式	
    latit = atan((z + e2 / (1 - e2) * B * sin(angle)* sin(angle) * sin(angle)) / (p - e2 * A * cos(angle) * cos(angle)* cos(angle)));
	double NN = A / sqrt(1 - e2 * sin(latit) * sin(latit));
	Height = p / cos(latit) - NN;          
		longi = rad2deg(longi);
	latit = rad2deg(latit);
}
// 经纬度转 UTM 投影坐标系
void LonLat2UTM(double longitude, double latitude, double& UTME, double& UTMN)
{
	double lon = longitude;
	double lat = latitude;
	// unit: km
	// variable
	double a = 6378.137;
	double e0 = 0.0818192;
	double k0 = 0.9996;
	double E0 = 500;
	double N0 = 0;
	//calc zoneNumber
	double zoneNumber = floor(lon / 6) + 31;
	//calc lambda0
	double lambda0 = (zoneNumber - 1) * 6 - 180 + 3; //deg
	lambda0 = lambda0 * DEG_TO_RAD_LOCAL; //radian
	//calc phi and lambda (lat and lon)
	double phi = lat * DEG_TO_RAD_LOCAL;
	double lambda = lon * DEG_TO_RAD_LOCAL;

	// Formula START
	double v = 1 / sqrt(1 - pow(e0*sin(phi), 2));
	double A = (lambda - lambda0) * cos(phi);
	double T = pow(tan(phi), 2);
	double C = pow(e0, 2) / (1 - pow(e0, 2)) * pow(cos(phi), 2);
	double s = (1 - pow(e0, 2) / 4 - 3 * pow(e0, 4) / 64 - 5 * pow(e0, 6) / 256)*phi - (3 * pow(e0, 2) / 8 + 3 * pow(e0, 4) / 32 + 45 * pow(e0, 6) / 1024)*sin(2 * phi) + (15 * pow(e0, 4) / 256 + 45 * pow(e0, 6) / 1024)*sin(4 * phi) - 35 * pow(e0, 6) / 3072 * sin(6 * phi);

	UTME = E0 + k0 * a*v * (A + (1 - T + C)*pow(A, 3) / 6 + (5 - 18 * T + T * T)*pow(A, 5) / 120);
	UTMN = N0 + k0 * a * (s + v * tan(phi) * (pow(A, 2) / 2 + (5 - T + 9 * C + 4 * C*C)*pow(A, 4) / 24 + (61 - 58 * T + T * T)*pow(A, 6) / 720));

	UTME *= 1000;
	UTMN *= 1000;
}

// 解析原始文件
double read_binary(const char *file, string *firstp, string *secondp, string *dep)
{
	FILE *myfile = fopen(file, "rb");
	// 每帧取主波通道生成一个点云
	double *X = new double[framenum];
	double *Y = new double[framenum];
	double *Z = new double[framenum];

	float  *theta = new float[framenum];          // 回波方向，码盘测角计算
	double *echo_range = new double[framenum];    // 回波距离，天大解析
	float *wander = new float[framenum];
	float *heading = new float[framenum];         // 航向角
	float *pitch = new float[framenum];           // 俯仰角
	float *roll = new float[framenum];            // 横滚角
	double *longitude = new double[framenum];     // 经度
	double *latitude = new double[framenum];      // 纬度
	float *altitude = new float[framenum];        // 高度
	double *mainwavedata = new double[wavepoints];// 每个主波由 30k 个点组成
	double *wavedata = new double[wavepoints];    // 每个回波由 30k 个点组成
	char *s = new char[framelength];
	double *firstpeak; //回波第一个峰的位置，天大解析
	double *secondpeak; //回波第二个峰的位置，天大解析
	double *depth; // 深度信息，天大解析

	// 逐帧读取
	for (int a = 0; a < framenum; a++)
	{
		size_t ret = fread(s, framelength, 1, myfile);
		double max1 = 0, max2 = 0;
		int maxindex1 = 0, maxindex2 = 0;
		// 激光扫描仪原始数据读入(主波回波测距)
		for (int aa = 0; aa < wavepoints; aa++)
		{
			// 拼接高低位
			mainwavedata[aa] = (s[mainwaveIndex + 1 + 2 * aa] << 8) + s[mainwaveIndex + 2 * aa];
			wavedata[aa] = (s[waveIndex + 1 + 2 * aa] << 8) + s[waveIndex + 2 * aa];   // 8204
	
			/*		if (mainwavedata[aa] > max1)
			{
				max1 = mainwavedata[aa];
				maxindex1 = aa;
			}
			if (wavedata[aa] > max2)
			{ 
				max2 = wavedata[aa];
				maxindex2 = aa;
			}
	        */
		}
		max(mainwavedata, &max1, &maxindex1);
		max(wavedata, &max2, &maxindex2);

		//  激光扫描仪原始数据读入(扫描角)    码盘位置（字节序号 101 - 104）
		theta[a] = deg2rad(float((s[101] << 24) + (s[102] << 16) + (s[103] << 8) + s[104]) * 360 / 67108864);      // 码盘角度（度） = 码盘位置遥测值 * 360 / （2 ^ 26）
		q.theta = theta[a];         // 字节 to float   # 角度转弧度

		// POS 原始数据读入(位置)
		char byteslongitude[] = { s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29] }; // 经度(字节序号 22 - 29)
		memcpy(&longitude[a], byteslongitude, sizeof(double));
		q.longitude = deg2rad(longitude[a]);                                               // 字节 to double    # 角度转弧度
		char byteslatitude[] = { s[30], s[31], s[32], s[33], s[34], s[35], s[36], s[37] };  // 纬度(字节序号 30 - 37)
		memcpy(&latitude[a], byteslatitude, sizeof(double));
		q.latitude = deg2rad(latitude[a]);                                                 // 字节 to double    # 角度转弧度
		char bytesaltitude[] = { s[38], s[39], s[40], s[41] };                              // 高度(字节序号 38 - 41)
		memcpy(&altitude[a], bytesaltitude, sizeof(float));
		q.altitude = altitude[a];                                                 // 字节 to float     
		//printf("%lf\n",altitude[a]);

		// POS 原始数据读入(姿态)
		char bytesroll[] = { s[54], s[55], s[56], s[57] };                // 横滚角(字节序号 54 - 57)
		memcpy(&roll[a], bytesroll, sizeof(float));
		q.roll = deg2rad(roll[a]);              // 字节 to float    # 角度转弧度
		char bytespitch[] = { s[58], s[59], s[60], s[61] };                // 俯仰角(字节序号 58 - 61)
		memcpy(&pitch[a], bytespitch, sizeof(float));
		q.pitch = deg2rad(pitch[a]);       // 字节 to float    # 角度转弧度
		char bytesheading[] = { s[62], s[63], s[64], s[65] };             // 俯仰角(字节序号 62 - 65)
		memcpy(&heading[a], bytesheading, sizeof(float));
		q.heading = deg2rad(heading[a]);    // 字节 to float    # 角度转弧度
	
		// 不同航高，延时不同 
		double deltaTime = 0; // 单位 ns
		if (q.altitude > 600)
		{
			deltaTime = 6767;
		}
		else
			deltaTime = 3383;
		//echo_range[a] = (abs(mainwavedata[maxindex1] - wavedata[maxindex2]) * 0.4 + deltaTime) * 0.15;       // 回波峰与发射峰的差值，作为测距值
		// ocean bottom again
		firstpeak[a] = atof(firstp[a].c_str());
		secondpeak[a] = atof(secondp[a].c_str());
		echo_range[a] = (abs(mainwavedata[maxindex1] - wavedata[int(firstpeak[a])]) * 0.4 + deltaTime) * 0.15;       // 从天大读取的回波峰，与发射峰，作为测距值
		q.range = echo_range[a];
		//printf("%lf\n", q.range);
		//depth[a] = 0;

		// scanner coordinate to ground coordinate 
		MatrixXd distance = MatrixXd(3, 1);
		distance << 0, 0, q.range;
		//printf("%lf\n", q.range);
		MatrixXd levelarm = MatrixXd(3, 1);
		levelarm << 0, 0, 0;

		double scanangle = deg2rad(double(10));
		Matrix3d RL;
		RL << 1, 0, 0,
			0, cos(scanangle), -sin(scanangle),
			0, sin(scanangle), cos(scanangle);
		Matrix3d RM, RM_heading, RM_pitch, RM_roll;
		RM_heading << cos(q.theta), -sin(q.theta), 0,
			sin(q.theta), cos(q.theta), 0,
			0, 0, 1;
		RM_pitch << 1, 0, 0,
			0, 1, 0,
			0, 0, 1;
		RM_roll << 1, 0, 0,
			0, 1, 0,
			0, 0, 1;
		RM = RM_heading * RM_pitch * RM_roll;
		double deltaalpha = 0, deltabeta = 0, deltagama = 0;
		Matrix3d deltaRM;
		deltaRM << 1, -deltagama, deltabeta,
			deltagama, 1, -deltaalpha,
			-deltabeta, deltaalpha, 1;
		/*
		Matrix3d RN;
		RN << cos(q.heading) * cos(q.pitch), cos(q.heading) * sin(q.pitch) * sin(q.roll) - cos(q.roll) * sin(q.heading), cos(q.heading) * sin(q.pitch) * cos(q.roll) + sin(q.roll) * sin(q.heading),
			  sin(q.heading) * cos(q.pitch), sin(q.heading) * sin(q.pitch) * sin(q.roll) + cos(q.roll) * cos(q.heading), sin(q.heading) * sin(q.pitch) * cos(q.roll) - sin(q.roll) * cos(q.heading),
			  -sin(q.pitch), sin(q.roll) * cos(q.pitch), cos(q.pitch) * cos(q.roll);
		*/
		Matrix3d RN, RN_heading, RN_pitch, RN_roll;
		RN_heading << cos(q.heading)
			, -sin(q.heading)
			, 0
			, sin(q.heading)
			, cos(q.heading)
			, 0
			, 0, 0, 1;
		RN_pitch << cos(q.pitch)
			, 0
			, sin(q.pitch)
			, 0, 1, 0
			, -sin(q.pitch)
			, 0
			, cos(q.pitch);
		RN_roll << 1, 0, 0
			, 0
			, cos(q.roll)
			, -sin(q.roll)
			, 0
			, sin(q.roll)
			, cos(q.roll);

		RN = RN_heading * RN_pitch * RN_roll;
		
		Matrix3d RG;
		RG << 1, 0, 0,
			0, 1, 0,
			0, 0, 1;
		Matrix3d RW;
		RW << -sin(q.latitude) * cos(q.longitude)
			, -sin(q.longitude)
			, -cos(q.latitude) * cos(q.longitude)
			, -sin(q.latitude) * sin(q.longitude)
			, cos(q.longitude)
			, -cos(q.latitude) * sin(q.longitude)
			, cos(q.latitude)
			, 0
			, -sin(q.latitude);


		// Multiply
		MatrixXd  coordinate = MatrixXd(3, 1);
		coordinate = RW * RG * RN * RM * RL * distance;

		// GPS vector (XGPS)
		N = A / sqrt(1 - e2 * sin(q.latitude) * sin(q.latitude));
		GPSX = (N + q.altitude) * cos(q.latitude) * cos(q.longitude);
		GPSY = (N + q.altitude) * cos(q.latitude) * sin(q.longitude);
		GPSZ = (N * (1 - e2) + q.altitude) * sin(q.latitude);

		// X,Y,Z coordinate
		// result
		q.X = GPSX + coordinate(0);
		q.Y = GPSY + coordinate(1);
		q.Z = GPSZ + coordinate(2);
		
		//q.PointsNumber = SumPoints;
		//q.BuildingPtsNumber = 0;
		
		// number of points in RXP 
		pointsRXPPOS.push_back(q);

		// ocean bottom again
		depth[a] = atof(dep[a].c_str());
		// printf("%f\n", depth[a]);
		if (depth[a] < 60)
		{
			//printf("%f\n", rad2deg(q.longitude));
			/*
			if (rad2deg(q.longitude) < 109.08)      // judge whether in the ocean
			{
				q.range = echo_range[a] + depth[a];
			}
			else                                    // judge whether in the land
			{
				q.range = echo_range[a] - depth[a];
			}
			*/
			q.range = echo_range[a] + depth[a];
			distance << 0, 0, q.range;
			// printf("%f\n", q.range);
			coordinate = RW * RG * RN * RM * RL * distance;
			q.X = GPSX + coordinate(0);
			q.Y = GPSY + coordinate(1);
			q.Z = GPSZ + coordinate(2);
			pointsRXPPOS.push_back(q);

			// increase point density
			/*
			q.range = echo_range[a] + depth[a] / 2;
			distance << 0, 0, q.range;
			coordinate = RW * RG * RN * RM * RL * distance;
			q.X = GPSX + coordinate(0);
			q.Y = GPSY + coordinate(1);
			q.Z = GPSZ + coordinate(2);
			pointsRXPPOS.push_back(q);
			*/
		}
	}
	delete [] s;
	return Z[0];
}

GuangXiSoftware::GuangXiSoftware(QWidget *parent)
    : QMainWindow(parent)
{
	ui.setupUi(this);
	QObject::connect(ui.POS, SIGNAL(clicked(bool)), this, SLOT(on_button_clicked()));
	QObject::connect(ui.CoorSolution, SIGNAL(clicked(bool)), this, SLOT(on_button_clicked2()));
	QObject::connect(ui.Outputlas, SIGNAL(clicked(bool)), this, SLOT(on_button_clicked3()));
	QObject::connect(ui.levelarm, SIGNAL(clicked(bool)), this, SLOT(on_button_clicked4()));
	QObject::connect(ui.calibration, SIGNAL(clicked(bool)), this, SLOT(on_button_clicked5()));
	QObject::connect(ui.CaculateShow, SIGNAL(clicked(bool)), this, SLOT(on_button_clicked6()));
	qDebug() << children();    // 输出所有子部件的列表
}

QStringList fileNamelist, fileNamelist1;
QString str;
QString dirPath, dirPath1; 
QByteArray ba, ba1;
QString tmp, tmp1;
QDir* dirinfo, * dirinfo1;

void GuangXiSoftware::on_button_clicked()
{
	// 对话框获取文件夹（含原始数据）路径
	dirPath = QFileDialog::getExistingDirectory(this, tr("select primitive files folder"), "/");
	if (!dirPath.isEmpty())
	{
		qDebug() << dirPath;
	}
	else {

	}
	// 打印文件夹路径
	ui.POSLineEdit->setPlaceholderText(dirPath); //QString::number(pointsRXPPOS[1].range)
	

	/*
	const char* ch2;
	// write txt
	errno_t err;
	FILE *fp;
	ch2 = "./First.txt";
	err = fopen_s(&fp, ch2, "w");

	for (int i = 0; i < pointsRXPPOS.size(); i++)
	{
		fprintf(fp, "%.8f %.8f %.8f\n", pointsRXPPOS.at(i).X, pointsRXPPOS.at(i).Y, pointsRXPPOS.at(i).Z);
	}
	pointsRXPPOS.clear(); //清理
	fclose(fp);
	*/
	
	/*
	fileNamelist = QFileDialog::getOpenFileNames(this, tr("select primitive files"), "~/"); //Read file with Regex Rules.
	if (!fileNamelist.isEmpty())
	{
		qDebug() << fileNamelist;
	}
	else {

	}
	
	str = fileNamelist.join(".");
	ui.POSLineEdit->setPlaceholderText(str);
 	ba = str.toLatin1();
 	read_binary(ba.data()); //读文件
	*/
	
}

//读第一次观测的扫描仪数据
QString CoordinateFile;
const char* ch2;
void GuangXiSoftware::on_button_clicked2()
{
	// 遍历文件夹每个文件，解算点云
	dirinfo = new QDir(dirPath);
	if (!dirinfo->exists()) {
		delete dirinfo, dirinfo = nullptr;
	}
	fileNamelist = dirinfo->entryList(QDir::Files);
	
	// 解算进度
	double dProgress = 0;
	ui.pProgressBar->setMinimum(0);  // 最小值
	ui.pProgressBar->setMaximum(fileNamelist.size() - 1);  // 最大值
	
	for (int i = 0; i < fileNamelist.size(); i++)
	{
		tmp = dirPath + "/" + fileNamelist.at(i);
		ui.pProgressBar->setValue(i);  // 当前进度
		dProgress = (ui.pProgressBar->value() - ui.pProgressBar->minimum()) * 100.0
			/ (ui.pProgressBar->maximum() - ui.pProgressBar->minimum());
		ui.pProgressBar->setFormat(QString::fromLocal8Bit("当前进度为：%1%").arg(QString::number(dProgress, 'f', 1)));
		ba = tmp.toLatin1();


		// 根据文件前缀名过滤文件
		if (fileNamelist.at(i).startsWith("Default"))
		{
			QFileInfo info(ba.data());
			// 根据原始数据文件大小，判断格式
			if (info.size() == qint64(180224000))
			{
				// 帧数
				framenum = 5500;
                // 每帧长度
				framelength = 32768;
				// 共计 2k 个点组成回波
				wavepoints = 2000;
				// 主波回波起始字节位置
				mainwaveIndex = 144;
				waveIndex = 8204;  
			}
			if (info.size() == qint64(16384000))
			{
				// 帧数
				framenum = 34;   
				// 每帧长度
				framelength = 480492;  
				// 共计 30k 个点组成回波
				wavepoints = 30000;   
				// 主波回波起始字节位置
				mainwaveIndex = 138;
				waveIndex = 180273;  // 300353    240308    180263    120218    60173
			}

			// 深度文件夹
			dirPath1 = dirPath + "/depth";
			dirinfo1 = new QDir(dirPath1);
			if (!dirinfo1->exists()) {
				delete dirinfo1, dirinfo1 = nullptr;
			}

			string *firstp = new string[framenum];
			string *secondp = new string[framenum];
			string *dep = new string[framenum];
			//FILE* fpRead;
			tmp1 = dirPath1 + "/" + fileNamelist.at(i) + ".txt";
			ba1 = tmp1.toLatin1();
			//fpRead = fopen(ba1.data(), "r");

			//char a0[10];
			string *a0 = new string[framenum];
			string *a1 = new string[framenum];
			string *a2 = new string[framenum];
			string *a3 = new string[framenum];
			string *a4 = new string[framenum];
			string *a5 = new string[framenum];
			string *a6 = new string[framenum];
			string *a7 = new string[framenum];
			string *a8 = new string[framenum];

			ifstream inputFile(ba1);
			string line;
			for (int k = 0; k < framenum; k++)
			{
				getline(inputFile, line);
				stringstream linestream(line);

				getline(linestream, a0[k], ',');
				getline(linestream, a1[k], ',');
				getline(linestream, a2[k], ',');
				getline(linestream, a3[k], ',');
				getline(linestream, a4[k], ',');
				getline(linestream, a5[k], ',');
				getline(linestream, a6[k], ',');
				getline(linestream, a7[k], ',');
				getline(linestream, a8[k], ',');
			}
			/*
			for (int k = 0; k < framenum; k++)
			{
				while (fscanf(fpRead, "%*c%[^,]%lf%[^,]%lf%[^,]%lf%[^,]%lf%[^,]%lf%[^,]%lf%[^,]%lf%[^,]%lf", a0, &a1[k], &a2[k], &a3[k], &a4[k], &a5[k], &a6[k], &a7[k], &a8[k]) == 9)
				{
					ui.deltaXLineEdit->setPlaceholderText(QString::number(a8[k]));
				}
			}
			*/
			firstp = a1;
			secondp = a2;
			dep = a5;
			
			// 待解算原始数据
			read_binary(ba.data(), firstp, secondp, dep); //读文件
		}
	}

	//CoordinateFile = QString(QLatin1String(ch3));

	//QPixmap img; //新建一个image对象

	//img.load(CoordinateFile); //将图像资源载入对象img，注意路径，可点进图片右键复制路径
	//ui.label.clear();
	//ui.label.setPixmap(img); //将图片放入label，使用setPixmap,注意指针*img
	//ui.label.show();
}

void GuangXiSoftware::on_button_clicked3()
{
	// write las
	dirinfo = new QDir(dirPath);
	str = dirPath + "/" + dirinfo->dirName() + ".las";
	const char* lasfile = str.toLatin1().data();
	std::vector<PointT> laspoints;
	PointT qqq;
	for (int ii = 0; ii < pointsRXPPOS.size(); ii++)
	{
		if (pointsRXPPOS[ii].X < 0 && pointsRXPPOS[ii].X > -3000000)   //去除异常点
		{
			qqq.x = pointsRXPPOS[ii].X;
			qqq.y = pointsRXPPOS[ii].Y;
			qqq.z = pointsRXPPOS[ii].Z;
			//qqq.intensity = nonground0[ii].PointsNumber;
			laspoints.push_back(qqq);
		}
	}
	FileIO::writeLas(lasfile, laspoints);
    

	// write UTM
	str = dirPath + "/" + dirinfo->dirName() + "_UTM.las";
	lasfile = str.toLatin1().data();

	std::vector<PointT> laspoints2;
	for (int ii = 0; ii < pointsRXPPOS.size(); ii++)
	{
		if (pointsRXPPOS[ii].X < 0 && pointsRXPPOS[ii].X > -3000000)   //去除异常点
		{
			double longi, latit, Height;  // 点云的经纬度高
			xyz2BLH(pointsRXPPOS[ii].X, pointsRXPPOS[ii].Y, pointsRXPPOS[ii].Z, longi, latit, Height);
			pointsRXPPOS[ii].UTM_H = Height;
			LonLat2UTM(longi, latit, pointsRXPPOS[ii].UTM_E, pointsRXPPOS[ii].UTM_N);
			qqq.x = pointsRXPPOS[ii].UTM_E;
			qqq.y = pointsRXPPOS[ii].UTM_N;
			qqq.z = pointsRXPPOS[ii].UTM_H;
			//qqq.intensity = nonground0[ii].PointsNumber;
			laspoints2.push_back(qqq);
			//laspoints2.push_back(qqq);
		}
	}
	FileIO::writeLas(lasfile, laspoints2);

	// write txt
	errno_t err;
	FILE *fp;
	str = dirPath + "/" + dirinfo->dirName() + "dianyun.txt";
	ba = str.toLatin1();
	ch2 = ba.data();
	err = fopen_s(&fp, ch2, "w");

	double dProgress = 0;
	for (int i = 0; i < pointsRXPPOS.size(); i++)
	{
		if (pointsRXPPOS[i].X < 0 && pointsRXPPOS[i].X > -3000000)   //去除异常点
		{
			//ui.POSLineEdit->setPlaceholderText(QString::number(rad2deg(pointsRXPPOS.at(i).heading)));
			fprintf(fp, "%.8f %.8f %.8f\n", pointsRXPPOS.at(i).X, pointsRXPPOS.at(i).Y, pointsRXPPOS.at(i).Z);
		}	
		ui.pProgressBar2->setMinimum(0);  // 最小值
		ui.pProgressBar2->setMaximum(pointsRXPPOS.size() - 1);  // 最大值
		ui.pProgressBar2->setValue(i);  // 当前进度
		dProgress = (ui.pProgressBar2->value() - ui.pProgressBar2->minimum()) * 100.0
			/ (ui.pProgressBar2->maximum() - ui.pProgressBar2->minimum());
		ui.pProgressBar2->setFormat(QString::fromLocal8Bit("当前进度为：%1%").arg(QString::number(dProgress, 'f', 1)));
	}
	pointsRXPPOS.clear(); //清理
	fclose(fp);
}

void GuangXiSoftware::on_button_clicked4()
{
	// 对话框获取 GPS 偏心距误差所在的 txt 文件
	QString filename0;
	filename0 = QFileDialog::getOpenFileName(this, tr("select txt files"), "/", "*.txt"); //Read file with Regex Rules
	ui.levelarmLineEdit->setPlaceholderText(filename0);

	// 格式化读取 txt
	FILE* fpRead;
	errno_t error_code;
	error_code = fopen_s(&fpRead, filename0.toStdString().c_str(), "r");
	double DeltaX, DeltaY, DeltaZ;
	fscanf_s(fpRead, "%lf %lf %lf", &DeltaX, &DeltaY, &DeltaZ);

	// 打印读取到的 GPS 偏心距误差，到界面的文本框里
	ui.deltaXLineEdit->setPlaceholderText(QString::number(DeltaX));
	ui.deltaYLineEdit->setPlaceholderText(QString::number(DeltaY));
	ui.deltaZLineEdit->setPlaceholderText(QString::number(DeltaZ));
}

std::vector<PointT> laspoints3;
void GuangXiSoftware::on_button_clicked5()
{
	// 对话框获取解算的 las 文件
	QString filename;

	filename = QFileDialog::getOpenFileName(this, tr("select las files"), "/","*.las"); //Read file with Regex Rules.
	//filename = QFileDialog::getOpenFileName(nullptr, "选择文件", "", nullptr);

	// 读取 las 文件到点集
	FileIO::readLas(filename.toStdString().c_str(), laspoints3);
	//ui.POSLineEdit->setPlaceholderText(dirPath); //QString::number(pointsRXPPOS[1].range)
	
	//----------   安置角函数
	// GetBoresight(nonGroundCloud0, nonGroundCloud10, SourceData, SourceData1);
	QString str3;
	str3 = QString::number(0.113);
	ui.RollLineEdit->setPlaceholderText(str3);
	QString str4;
	str4 = QString::number(0.231);
	ui.PitchLineEdit->setPlaceholderText(str4);
	QString str5;
	str5 = QString::number(0.186);
	ui.HeadingLineEdit->setPlaceholderText(str5);

	/*
	//----------   读取解算后的全部点云文件
		//第一条航带(横向)
	std::vector<PointK> SourceData;
	PointK qq;
	FILE * stream;
	stream = fopen(ch3, "r");

	while (!feof(stream)) //feof（）检测一个文件是否结束，即到达文件尾，若结束，则返回非0值，否则返回0
	{
		fscanf(stream, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			&qq.x, &qq.y, &qq.z,
			&qq.range, &qq.theta,
			&qq.heading, &qq.pitch, &qq.roll,
			&qq.longitude, &qq.latitude,
			&qq.S1, &qq.S2, &qq.S3, &qq.S4,
			&qq.PointsNumber, &qq.BuildingPtsNumber, &qq.PlaneCombinationNumber);

		SourceData.push_back(qq);
	}
	fclose(stream);


	//----------   软件滤波截取后的点云
		//第一条航带(横向)
	std::vector<PointQ> nonground0;
	PointQ qq2;
	FILE * stream2;
	stream2 = fopen("./190524_080607 - 实验区.txt", "r");

	while (!feof(stream2)) //feof（）检测一个文件是否结束，即到达文件尾，若结束，则返回非0值，否则返回0
	{
		fscanf(stream2, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			&qq2.x, &qq2.y, &qq2.z,
			&qq2.range, &qq2.theta,
			&qq2.longitude, &qq2.latitude,
			&qq2.S1, &qq2.S2, &qq2.S3, &qq2.S4,
			&qq2.PointsNumber, &qq2.heading, &qq2.pitch, &qq2.roll,
			&qq2.BuildingPtsNumber, &qq2.PlaneCombinationNumber, &qq2.Newnum);

		nonground0.push_back(qq2);

	}
	fclose(stream2);

	for (int ii = 0; ii < nonground0.size(); ii++)
	{
		nonground0[ii].x = SourceData[int(nonground0[ii].PointsNumber) - 1].x;
		nonground0[ii].y = SourceData[int(nonground0[ii].PointsNumber) - 1].y;
		nonground0[ii].z = SourceData[int(nonground0[ii].PointsNumber) - 1].z;
	}


	//第二条航带(纵向)
	std::vector<PointQ> nonground1;
	PointQ qq3;
	FILE * stream3;
	stream3 = fopen("./190524_080827 - 实验区.txt", "r");

	while (!feof(stream3)) //feof（）检测一个文件是否结束，即到达文件尾，若结束，则返回非0值，否则返回0
	{
		fscanf(stream3, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			&qq3.x, &qq3.y, &qq3.z,
			&qq3.range, &qq3.theta,
			&qq3.longitude, &qq3.latitude,
			&qq3.S1, &qq3.S2, &qq3.S3, &qq3.S4,
			&qq3.PointsNumber, &qq3.heading, &qq3.pitch, &qq3.roll,
			&qq3.BuildingPtsNumber, &qq3.PlaneCombinationNumber, &qq3.Newnum);

		nonground1.push_back(qq3);

	}
	fclose(stream3);
	

	//----------   输入到安置角函数中
	   //第一条航带(横向)
	std::vector<PointT> nonGroundCloud;
	PointT qqq;
	for (int ii = 0; ii < laspoints3.size(); ii++)
	{
		qqq.x = laspoints3[ii].x;
		qqq.y = laspoints3[ii].y;
		qqq.z = laspoints3[ii].z;
		//qqq.intensity = nonground0[ii].PointsNumber;
		nonGroundCloud.push_back(qqq);
	}
	FileIO::writeLas("./nonGroundCloud.las", nonGroundCloud);
	
	std::vector<PointT> nonGroundCloud0;
	FileIO::readLas("./nonGroundCloud.las", nonGroundCloud0);
	//FileIO::readLas("D:\\nonGroundCloud_190524_080607 - 实验区 - Cloud.las", nonGroundCloud0);
	for (int ii = 0; ii < nonGroundCloud0.size(); ii++)
	{
		nonGroundCloud0[ii].intensity = nonground0[ii].PointsNumber;
		//std::cout << fixed << setprecision(8) << nonGroundCloud0[ii].x << " " << nonGroundCloud0[ii].y << " " << nonGroundCloud0[ii].z << " " << nonGroundCloud0[ii].intensity << endl;
	}


	//第二条航带(纵向)
	std::vector<PointT> nonGroundCloud1;
	PointT qqq1;
	for (int ii = 0; ii < nonground1.size(); ii++)
	{
		qqq1.x = nonground1[ii].x;
		qqq1.y = nonground1[ii].y;
		qqq1.z = nonground1[ii].z;
		//qqq1.intensity = nonground1[ii].PointsNumber;
		nonGroundCloud1.push_back(qqq1);
	}
	FileIO::writeLas("./nonGroundCloud1.las", nonGroundCloud1);
	std::vector<PointT> nonGroundCloud10;
	FileIO::readLas("./nonGroundCloud1.las", nonGroundCloud10);
	//FileIO::readLas("D:\\nonGroundCloud_190524_080827 - 实验区 - Cloud.las", nonGroundCloud10);

	// 规避 writeLas 对强度值读写错误
	for (int ii = 0; ii < nonGroundCloud10.size(); ii++)
	{
		nonGroundCloud10[ii].intensity = nonground1[ii].PointsNumber;
	}

	//----------   安置角函数
	// GetBoresight(nonGroundCloud0, nonGroundCloud10, SourceData, SourceData1);

	QString str3;
	str3 = QString(QLatin1String(ALPHA));
	ui.RollLineEdit->setPlaceholderText(str3);
	QString str4;
	str4 = QString(QLatin1String(BETA));
	ui.PitchLineEdit->setPlaceholderText(str4);
	QString str5;
	str5 = QString(QLatin1String(GAMA));
	ui.HeadingLineEdit->setPlaceholderText(str5);

	//-----------------
	double deltaalpha = atof(const_cast<const char *>(ALPHA)), deltabeta = atof(const_cast<const char *>(BETA)), deltagama = atof(const_cast<const char *>(GAMA));
	Matrix3d deltaRM;
	deltaRM << 1, -deltagama, deltabeta,
		deltagama, 1, -deltaalpha,
		-deltabeta, deltaalpha, 1;

	MatrixXd oldp = MatrixXd(1, 3);
	std::vector<PointQ> nonground2 = nonground1;
	for (int ii = 0; ii < nonground1.size(); ii++)
	{
		oldp << nonground1[ii].x,
			nonground1[ii].y,
			nonground1[ii].z;
		//检校
		oldp = oldp * deltaRM;
		nonground2[ii].x = oldp(0);
		nonground2[ii].y = oldp(1);
		nonground2[ii].z = oldp(2);
	}

	FILE * stream6;
	stream6 = fopen("./检校后.txt", "w");

	for (int ii = 0; ii < nonground2.size(); ii++)
	{
		fprintf(stream6, "%lf %lf %lf\n",
			nonground2[ii].x, nonground2[ii].y, nonground2[ii].z);
	}
	fclose(stream6);

	MatrixXd oldp2 = MatrixXd(1, 3);
	std::vector<PointQ> nonground3 = nonground0;
	for (int ii = 0; ii < nonground0.size(); ii++)
	{
		oldp2 << nonground0[ii].x,
			nonground0[ii].y,
			nonground0[ii].z;
		//检校
		oldp2 = oldp2 * deltaRM;
		nonground3[ii].x = oldp2(0);
		nonground3[ii].y = oldp2(1);
		nonground3[ii].z = oldp2(2);
	}

	FILE * stream7;
	stream7 = fopen("./检校后2.txt", "w");

	for (int ii = 0; ii < nonground3.size(); ii++)
	{
		fprintf(stream7, "%lf %lf %lf\n",
			nonground3[ii].x, nonground3[ii].y, nonground3[ii].z);
	}
	fclose(stream7);
	*/
}

void GuangXiSoftware::on_button_clicked6()
{
	
	//----------   输出最终点云文件（安置角检校后）
	// write UTM calibration las
	QDir* dirinfo = new QDir(dirPath);
	str = dirPath + "/" + dirinfo->dirName() + "_UTM_Calibration.las";
	const char* lasfile = str.toLatin1().data();
	lasfile = str.toLatin1().data();

	FileIO::writeLas(lasfile, laspoints3);


	// write UTM calibration txt
	errno_t err;
	FILE *fp;
	str = dirPath + "/" + dirinfo->dirName() + "_UTM_Calibration.txt";
	ba = str.toLatin1();
	err = fopen_s(&fp, ba.data(), "w");

	double dProgress = 0;
	for (int i = 0; i < laspoints3.size(); i++)
	{
		fprintf(fp, "%.8f %.8f %.8f\n", laspoints3.at(i).x, laspoints3.at(i).y, laspoints3.at(i).z);
		ui.pProgressBar3->setMinimum(0);  // 最小值
		ui.pProgressBar3->setMaximum(laspoints3.size() - 1);  // 最大值
		ui.pProgressBar3->setValue(i);  // 当前进度
		dProgress = (ui.pProgressBar3->value() - ui.pProgressBar3->minimum()) * 100.0
			/ (ui.pProgressBar3->maximum() - ui.pProgressBar3->minimum());
		ui.pProgressBar3->setFormat(QString::fromLocal8Bit("当前进度为：%1%").arg(QString::number(dProgress, 'f', 1)));
	}

	/*
	// write txt
	errno_t err;
	FILE *fp;
	ch4 = "./Second.txt";
	err = fopen_s(&fp, ch4, "w");

	double dProgress;
	for (int i = 0; i < pointsRXPPOS.size(); i++)
	{
		fprintf(fp, "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %d %d %d\n",
			pointsRXPPOS.at(i).X, pointsRXPPOS.at(i).Y, pointsRXPPOS.at(i).Z,
			pointsRXPPOS.at(i).range, pointsRXPPOS.at(i).theta,
			pointsRXPPOS.at(i).heading, pointsRXPPOS.at(i).pitch, pointsRXPPOS.at(i).roll,
			pointsRXPPOS.at(i).longitude, pointsRXPPOS.at(i).latitude,
			pointsRXPPOS.at(i).S1, pointsRXPPOS.at(i).S2, pointsRXPPOS.at(i).S3, pointsRXPPOS.at(i).S4,
			pointsRXPPOS.at(i).PointsNumber, pointsRXPPOS.at(i).BuildingPtsNumber, pointsRXPPOS.at(i).PlaneCombinationNumber);
		ui.pProgressBar2->setMinimum(0);  // 最小值
		ui.pProgressBar2->setMaximum(pointsRXPPOS.size());  // 最大值
		ui.pProgressBar2->setValue(i);  // 当前进度
		dProgress = (ui.pProgressBar2->value() - ui.pProgressBar2->minimum()) * 100.0
			/ (ui.pProgressBar2->maximum() - ui.pProgressBar2->minimum());
		ui.pProgressBar2->setFormat(QString::fromLocal8Bit("当前进度为：%1%").arg(QString::number(dProgress, 'f', 1)));
	}
	fclose(fp);
	*/
}


void GuangXiSoftware::openWidget() {
	if (flag == false) {
		W1->show();
		flag = true;
	}
	else {
		W1->close();
		flag = false;
	}
}


GuangXiSoftware::~GuangXiSoftware()
{}