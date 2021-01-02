#ifndef _COMP_H_
#define _COMP_H_

#include <armadillo>
#include <vector>
#include <assert.h>
#include <iomanip>


namespace civ {
    typedef enum {
        k,
        m
    } global_matrix;
    class element {
    public:
    
        element(float E, float I, float A, float L, float rho,
            float x_1, float y_1, float x_2, float y_2, float theta);
        ~element(); 

        std::array<float, 2> x1y1();
        std::array<float, 2> x2y2();
        std::array<std::array<float, 2>, 2> pts();
        
        // TO DO: DEPRECATE THESE AND CHANGE THIS TO A SIMPLE ARRRAY OR ST 
        float x1() { return x1_; }
        float y1() { return y1_; }
        float x2() { return x2_; }
        float y2() { return y2_; }
        int x1_g() { return x1_g_; }
        int y1_g() { return y1_g_; }
        int r1_g() { return r1_g_; }
        int x2_g() { return x2_g_; }
        int y2_g() { return y2_g_; }
        int r2_g() { return r2_g_; }

        void setGlobalxyr(int x, int y, int r, int pside); 
        int getMaxElementGlobalDim();
        arma::mat& ConvertToGlobalSize(global_matrix k, int g_size); 

        // returns 6x6 globally transformed (K = T^T * k * T) matrix
        arma::mat& getGlobal(global_matrix matrix);

        // returns g_size * g_size globally transformed matrix for frame addition
        arma::mat getGlobalFrame(global_matrix matrix, int global_size);
        arma::mat& k_g() { return global_k; }
        arma::mat& m_g() { return global_m; }
        

    private: 
        void local_k(); 
        void trans();
        void local_m();
        void local2global(); 
        arma::mat k_l, m_l, t_, global_k, global_m;  
        float E_, I_, A_, L_, rho_, x1_, y1_, x2_, y2_, theta_;

        int x1_g_ = 0, y1_g_ = 0, r1_g_ = 0, 
            x2_g_ = 0, y2_g_ = 0, r2_g_ = 0;  
    };

    class frame {
    public:
        frame(); 
        ~frame(); 
        void add_element(element* elem);
        void global_k(); 
        void printPointDataForStructure();
        void printPoints(std::array<float, 2> point);
        void printGCS(); 

    private:
        std::vector<element*> elements_; 

        void checkElementBounds(element* elem); 

        // Global Coordinate System 
        std::vector<std::array<float, 2>> FindCoincidentPoints(); 
        void NumberSharedPoints(int xg, int yg, int rg, std::array<float, 2> shared_pt); 
        void CreateGCS();
        int GetMaximumGCSDim(); 

        arma::mat k_g_;          // global k matrix 
        arma::mat m_g_;          // global m matrix

        std::array<int, 3> g_coord = { 1, 2, 3 };
        bool GCS_created = false; 
        friend class element; 
        
    };
}

#endif // !_COMP_H_