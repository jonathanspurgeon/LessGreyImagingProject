#ifndef PARTICLE_H
#define PARTICLE_H

class particle
{
public:
    particle();
    int getId() const;
    void setId(int value);

    int getPoreID() const;
    void setPoreID(int value);

    double getXCoordinate() const;
    void setXCoordinate(double value);

    double getYCoordinate() const;
    void setYCoordinate(double value);

    double getZCoordinate() const;
    void setZCoordinate(double value);

    double getPorePosition() const;
    void setPorePosition(double value);

    bool getClosed() const;
    void setClosed(bool value);

private:
    int id;
    int poreID;
    double xCoordinate;
    double yCoordinate;
    double zCoordinate;
    double porePosition;
    bool closed;
};

#endif // PARTICLE_H
