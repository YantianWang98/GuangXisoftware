#include "GuangXiSoftware.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    GuangXiSoftware w;
    w.setWindowTitle(QString::fromLocal8Bit("µ„‘∆Ω‚À„"));
	w.show();
	
    return a.exec();
}