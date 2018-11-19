#ifndef SKELETON_H
#define SKELETON_H
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#ifdef _WIN32
#include "GL/glut.h"
#else
#include <GL/freeglut.h>
#endif
#endif

#include "simpleMath.h"
#include "Quaternion.h"
#include <iostream>


struct Joint
{
    int parentID;
    Vec3 position;
    Vec2 screenCoord;
    bool isHovered;
    bool isPicked;
    
    float local_t[16];
	Quaternion localQuaternion;
    float global_t[16];
    
    Joint()
    {
        parentID = -1;
        isHovered = false;
        isPicked = false;
        
        for(int i=0; i<16; ++i)
            local_t[i]=0.0;
        local_t[0]=local_t[5]=local_t[10]=local_t[15]=1.0;
        
        for(int i=0; i<16; ++i)
            global_t[i]=0.0;
        global_t[0]=global_t[5]=global_t[10]=global_t[15]=1.0;
    }

	void loadFromLocalQuaternion(Quaternion newLocalQuaternion) {
		localQuaternion = newLocalQuaternion;
		float* newLocal_t = newLocalQuaternion.toFloatMatrix();
		for(int i =0; i < 16; i++)
		{
			local_t[i] = newLocal_t[i];
		}
	}

	void setLocalTransform(float newLocal_t[16]) {
		for (int i = 0; i < 16; ++i) {
			local_t[i] = newLocal_t[i];
		}
	}
};

class Skeleton
{
private:
    /*Update screen coordinates of joints*/
    void updateScreenCoord();
    
public:
    std::vector<Joint> joints;
    /*True if the skeleton has a joint selected*/
    bool hasJointSelected;   
    Skeleton(){hasJointSelected = false;};
    /*
     * Load Skeleton file
     */
    void loadSkeleton(std::string skelFileName);

    /*
     * Load animation file
     */
    void loadAnimation(std::string skelFileName);

    /*
     * Draw skeleton with OpenGL
     */
    void glDrawSkeleton();

    /*
     * Check if any joint is hovered by given mouse coordinate
     */
    void checkHoveringStatus(int x, int y);

    void release();
    
    void updateGlobal();
    
    void addRotation(float* q);
    
    void selectOrReleaseJoint();
};

#endif
