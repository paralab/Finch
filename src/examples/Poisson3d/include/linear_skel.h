    #ifndef DENDRO_5_0_LINEAR_SKEL_H
    #define DENDRO_5_0_LINEAR_SKEL_H

    #include "oda.h"
    #include "feVector.h"
    #include "Genfunction.h"

    namespace FinchDendroSkeleton
    {
        class RHSVec : public feVector<RHSVec>{

        private:

            double * imV1;
            double * imV2;
            // function for boundary
	        std::function<void(double,double,double,double*)> bdry_function;

        public:
            RHSVec(ot::DA* da,unsigned int dof=1);
            ~RHSVec();
            
            /**@biref elemental compute vec for rhs*/
            virtual void elementalComputVec(const VECType* in,VECType* out, double*coords=NULL,double scale=1.0);
            
            /**@brief set boundary function*/	
            void setBdryFunction(std::function<void(double,double,double,double*)> bdry);
            
            bool preComputeVec(const VECType* in,VECType* out, double scale=1.0);
            
            bool postComputeVec(const VECType* in,VECType* out, double scale=1.0);
            
            /**@brief octree grid x to domin x*/
            double gridX_to_X(double x);
            /**@brief octree grid y to domin y*/
            double gridY_to_Y(double y);
            /**@brief octree grid z to domin z*/
            double gridZ_to_Z(double z);
        };
    }

    #endif //DENDRO_5_0_LINEAR_SKEL_H
    