#ifndef WORKER_H
#define WORKER_H

#include <QObject>
#include <iostream>
#include "network.h"

class worker : public QObject
{
    Q_OBJECT
public:
    explicit worker(network* net, int j=0, QObject *parent = 0);

    int getJob() const;
    void setJob(int value);

public slots:
    void process();

signals:
    void finished();

private:
    network* n;
    int job;

};

#endif // WORKER_H
