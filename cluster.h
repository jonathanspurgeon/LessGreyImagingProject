#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <map>

class cluster
{
public:
    cluster(int);
    int getLabel() const;
    void setLabel(int);

    bool getInlet() const;
    void setInlet(bool value);

    bool getOutlet() const;
    void setOutlet(bool value);

    bool getSpanning() const;
    void setSpanning(bool value);

private:
    int id;
    int label;
    bool inlet;
    bool outlet;
    bool spanning;
};

#endif // CLUSTER_H
