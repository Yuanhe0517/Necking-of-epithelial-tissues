/*!
 * \file constraint_tensile.hpp
 * \author Yuan He, heyuan@westlake.edu.cn
 * \date 4-Jun-2024
 * \brief Tensile constraint that sets force in x direction to zero
 */

#ifndef __CONSTRAINT_TENSILE_HPP__
#define __CONSTRAINT_TENSILE_HPP__

#include "constraint.hpp"

namespace VMTutorial
{
    class ConstraintTensile : public Constraint
    {
    public:
        ConstraintTensile() {}
        
        Vec apply(const Vertex<Property> &v, const Vec &f) override
        {
            Vec new_force = f;
            //new_force.x = 0.0;  // 设置x方向上的力为零
           	//cout<<new_force<<endl;
            return new_force;
        }

        Vec apply(const Vec &f) override
        {
            Vec new_force = f;
            new_force.x = 0.0;  // 设置x方向上的力为零
            return new_force;
        }
    };
}
#endif