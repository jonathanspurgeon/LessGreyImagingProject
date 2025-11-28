#include <QApplication>
#include <QThread>
#include <QObject>
#include "network.h"
#include "mainwindow.h"
#include "worker.h"

int main(int argc , char *argv[])
{

   QApplication a(argc,argv);

   MainWindow window;
   window.show();

   return a.exec();
}
