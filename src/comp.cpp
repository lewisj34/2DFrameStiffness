
#include "comp.h"
#include "macros.h"

using namespace civ;

element::element(float E, float I, float A, float L, float rho,
    float x1, float y1, float x2, float y2, float theta) 
    : E_(E), I_(I), A_(A), L_(L), rho_(rho),
    x1_(x1), y1_(y1), x2_(x2), y2_(y2), theta_(theta) {
    if (E <= 0 || I <= 0 || A <=0 || L <= 0 || rho <= 0) {
        printf("E, I, A, L, and rho must all be greater than 0\n");
        std::exit(EXIT_FAILURE); 
    }
#if DEBUG_ELEMENT
    std::cout << "*********************************************************" << std::endl;
    std::cout << "*****************       ELEMENT " << __COUNTER__ << "        ****************" << std::endl;
    std::cout << "*********************************************************" << std::endl;
#endif
    local_k(); 
    trans();
    local_m();
    local2global(); 
}

element::~element() {
    // do nothing... for now 
}

std::array<float, 2> element::x1y1() {
    std::array<float, 2> x1y1_out; 
    x1y1_out[0] = x1_;
    x1y1_out[1] = y1_;
    return x1y1_out; 
}

std::array<float, 2> element::x2y2() {
    std::array<float, 2> x2y2_out; 
    x2y2_out[0] = x2_;
    x2y2_out[1] = y2_;
    return x2y2_out; 
}

std::array<std::array<float, 2>, 2> element::pts() {
    std::array<std::array<float, 2>, 2> points;
    points[0] = x1y1();
    points[1] = x2y2();
    return points; 
}

void frame::printPoints(std::array<float, 2> point) { 
    for (int i = 0; i < point.size() - 1; i++)
        std::cout << "(" << point[i] << ", " << point[i+1] << ")";
}

void element::local_k() {

#if DEBUG_ELEMENT
    std::cout << "E = " << E_ << std::endl; 
    std::cout << "I = " << I_ << std::endl; 
    std::cout << "A = " << A_ << std::endl; 
    std::cout << "L = " << L_ << std::endl; 
    std::cout << "rho = " << rho_ << std::endl; 
#endif
    k_l.resize(6, 6); 
    k_l = { {A_ * pow(L_, 2) / I_, 0, 0, -A_ * pow(L_, 2) / I_, 0, 0},
            {0, 12, 6 * L_, 0, -12, 6 * L_},
            {0, 6 * L_, 4 * pow(L_,2), 0, -6*L_, 2 * pow(L_,2)},
            {-A_ * pow(L_, 2) / I_, 0, 0, A_ * pow(L_, 2) / I_, 0, 0},
            {0, -12, -6*L_, 0, 12, -6*L_},
            {0, 6* L_, 2 * pow(L_,2), 0, -6*L_, 4*pow(L_,2)} };

    k_l = k_l * (E_ * I_ / pow(L_,3));
#if DEBUG_ELEMENT
    std::cout << "k_l" << std::endl;
    std::cout << k_l << std::endl; 
    std::cout << std::endl;
#endif
}
void element::trans() {
    t_.resize(6, 6);

    theta_ *= pi / 180;
    float cos0_o;
    float sin0_o;
    float cos0 = std::cos(theta_);
    float sin0 = std::sin(theta_);

    bool modified_cos;
    bool modified_sin;
    if (cos0 < 1e-7) {
        cos0_o = cos0; 
        cos0 = 0;
        modified_cos = true;
    }
    if (sin0 < 1e-7) {
        sin0_o = sin0;
        sin0 = 0;
        modified_sin = true;
    }
#if DEBUG_ELEMENT
    std::cout << std::setprecision(4) << "cos0 = " << cos0 << std::endl;
    std::cout << "sin0 = " << sin0 << std::endl;
    if (modified_cos || modified_sin) {
        std::cout << "sin or cos modififed due to (sin0 || cos0) < 1e-7\n";
        std::cout << "\toriginal value of cos0 = " << cos0_o << std::endl;
        std::cout << "\toriginal value of sin0 = " << sin0_o << std::endl;
    }
#endif
    t_ = { {cos0, sin0, 0, 0, 0, 0},
        {-sin0, cos0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0},
        {0, 0, 0, cos0, sin0, 0},
        {0, 0, 0, -sin0, cos0, 0},
        {0, 0, 0, 0, 0, 1} };

#if DEBUG_ELEMENT
    std::cout << "t_" << std::endl;
    std::cout << t_ << std::endl; 
    std::cout << std::endl;
#endif
}


void element::local_m() {
    m_l.resize(6, 6); 
    m_l = { {2/6.f, 0, 0, 1/6.f, 0, 0},
        {0, 156/420.f, 22 * L_ /420.f, 0, 54 / 420.f, -13*L_/420.f},
        {0, 22 *L_/420.f, 4 * pow(L_,2)/420.f, 0, 13 *L_/420.f, -3*pow(L_,2)/420.f},
        {1/6.f, 0, 0, 2/6.f, 0, 0},
        {0, 54/420.f, 13*L_/420.f, 0, 156/420.f, -22*L_/420.f},
        {0, -13*L_/420.f, -3*pow(L_,2)/420.f, 0, -22*L_/420.f, 4 * pow(L_,2)/420.f} };
    m_l *= rho_ * A_ * L_; 
#if DEBUG_ELEMENT
    std::cout << "m_l" << std::endl;
    std::cout << m_l << std::endl; 
    std::cout << std::endl;
#endif
}

void element::local2global() {
    global_k = t_.t() * k_l * t_; 
	global_m = t_.t() * m_l * t_;

#if DEBUG_ELEMENT
    std::cout << "global_k" << std::endl;
    std::cout << global_k << std::endl; 
    std::cout << std::endl;
    std::cout << "global_m" << std::endl;
    std::cout << global_m << std::endl; 
    std::cout << std::endl;
#endif
}


frame::frame() {
    // do nothing... for now
}

frame::~frame() {
    if (elements_.size() > 0) {
        for (auto elem : elements_)
            delete elem; 
    }
}

void frame::add_element(element* elem) {
    checkElementBounds(elem); 
    elements_.push_back(elem); 
}


void frame::checkElementBounds(element* elem) {
    if (elements_.size() > 1) { 
        int num_pts = 0; // coincident points
        for (element* element : elements_) {
            if (elem->x1y1() == element->x1y1())
                num_pts += 1;
            if (elem->x1y1() == element->x2y2())
                num_pts += 1;
        }
        
        if (num_pts != 2) {
            std::cout << "Num coincident points = " << num_pts 
            << ". Only doing 2 coincident point cxns" 
            << std::endl; 
            std::exit(EXIT_FAILURE); 
        }
    }
    assert(elem->x1() >= 0 && elem->y1() >= 0);
    assert(elem->x2() >= 0 && elem->y2() >= 0);
}


void frame::printPointDataForStructure() {
    std::cout << "Printing point data for all elements in frame: " << std::endl; 
    std::cout << "i\t(x1,y1)\t\t(x2,y2)" << std::endl; 
    int i = 0; 
    for (element* elem : elements_) {
        std::cout << i << "\t(" << elem->x1() << ", " << elem->y1() << ")\t\t"  
        << "(" << elem->x2() << ", " << elem->y2() << ")" <<  std::endl; 
        i+=1;
    }
}

void frame::printGCS() {
    if (!GCS_created) {
        std::cout << "createGCS() not called, can't print GCS\n";
        std::exit(EXIT_FAILURE); 
    }
    std::cout << "Global coordinate system shown below" << std::endl; 
    std::cout << "i \tx1 \ty1 \tr1 \tx2 \ty2 \tr2" << std::endl; 
    for (int i = 0; i < elements_.size(); i++) {
        std::cout << i << "\t" << elements_[i]->x1_g()  
        << "\t" << elements_[i]->y1_g()  << "\t" << elements_[i]->r1_g() 
        << "\t" << elements_[i]->x2_g() << "\t" << elements_[i]->y2_g()  
        << "\t" << elements_[i]->r2_g() << std::endl; 
    }
}

// returns a vector of shared points found
std::vector<std::array<float, 2>> frame::FindCoincidentPoints() {
    // output container 
    std::vector<std::array<float, 2>> copts; 
    // iterate through all elements 
    for (int e = 0; e < elements_.size(); e++) {
        // now iterate through all elements other than current element (hence the e_c != e)
        for (int e_c = 0; e_c < elements_.size(); e_c++) {
            if (e != e_c) {
#if DEBUG_GCS
                std::cout << "Comparing element " << e << " with " << " element " << e_c << std::endl; 
#endif
                // point 1 -> point 1 
                if (elements_[e]->pts()[0] == elements_[e_c]->pts()[0]) {
#if DEBUG_GCS
                    std::cout << "\tE" << e << ", point 1 ";
                    printPoints(elements_[e]->pts()[0]); 
                    std::cout << " & E" << e_c << " point 1"; 
                    printPoints(elements_[e_c]->pts()[0]);
                    std::cout << " are the same." << std::endl; 
#endif
                    // shared point container (not necessary because of if statement)
                    std::array<float, 2> pt = elements_[e]->pts()[0];

                    // append coincident pt if it isn't already in the copoint vector
                    if (std::find(copts.begin(), copts.end(), pt) == copts.end())
                        copts.push_back(pt);
                }

                // point 1 -> point 2
                else if (elements_[e]->pts()[0] == elements_[e_c]->pts()[1]) {
#if DEBUG_GCS
                    std::cout << "\tE" << e << ", point 1 ";
                    printPoints(elements_[e]->pts()[0]); 
                    std::cout << " & E" << e_c << " point 2"; 
                    printPoints(elements_[e_c]->pts()[1]);
                    std::cout << " are the same." << std::endl; 
#endif
                    // shared point container (not necessary because of if statement)
                    std::array<float, 2> pt = elements_[e]->pts()[0];

                    // append coincident pt if it isn't already in the copoint vector
                    if (std::find(copts.begin(), copts.end(), pt) == copts.end())
                        copts.push_back(pt);
                }

                // point 2 -> point 1
                else if (elements_[e]->pts()[1] == elements_[e_c]->pts()[0]) {
#if DEBUG_GCS
                    std::cout << "\tE" << e << ", point 2 ";
                    printPoints(elements_[e]->pts()[1]); 
                    std::cout << " & E" << e_c << " point 1"; 
                    printPoints(elements_[e_c]->pts()[0]);
                    std::cout << " are the same." << std::endl;
#endif

                    // shared point container (not necessary because of if statement)
                    std::array<float, 2> pt = elements_[e]->pts()[1];

                    // append coincident pt if it isn't already in the copoint vector
                    if (std::find(copts.begin(), copts.end(), pt) == copts.end())
                        copts.push_back(pt);
                }

                // point 2 -> point 2
                else if (elements_[e]->pts()[1] == elements_[e_c]->pts()[1]) {
#if DEBUG_GCS
                    std::cout << "\tE" << e << ", point 2 ";
                    printPoints(elements_[e]->pts()[1]); 
                    std::cout << " & E" << e_c << " point 2"; 
                    printPoints(elements_[e_c]->pts()[1]);
                    std::cout << " are the same." << std::endl; 
#endif
                    // shared point container (not necessary because of if statement)
                    std::array<float, 2> pt = elements_[e]->pts()[1];
                    
                    // append coincident pt if it isn't already in the copoint vector
                    if (std::find(copts.begin(), copts.end(), pt) == copts.end())
                        copts.push_back(pt);
                }
                
            }
        }
    }
#if DEBUG_GCS
    std::cout << "Number of coincident points total: " 
    << copts.size() << "... Printing..." << std::endl;
    for (int i = 0; i < copts.size(); i++) {
        std::cout << "\t"; printPoints(copts[i]); std::cout << ", ";
    }
    std::cout << std::endl; 
#endif
    return copts;
}

void element::setGlobalxyr(int x, int y, int r, int pside) {
    assert(pside == 0 || pside == 1);
    if (pside == 0) { // x1y1r1
        x1_g_ = x;
        y1_g_ = y;
        r1_g_ = r;
    }
    else { // pside = 1 = x2y2r2
        x2_g_ = x;
        y2_g_ = y;
        r2_g_ = r;
    }
}

// finds maximum global dim amongst all globals: 
// x1_g_, y1_g_, r1_g_, x2_g_, y2_g_, r2_g_
int element::getMaxElementGlobalDim() {
    std::vector<int> gDims = { x1_g_, y1_g_, r1_g_, x2_g_, y2_g_, r2_g_ };
    int max_dim = *std::max_element(gDims.begin(), gDims.end());
    return max_dim;
}

arma::mat& element::getGlobal(global_matrix matrix) {
    if (matrix = global_matrix::k)
        return global_k; 
    else // matrix = mass matrix, m 
        return global_m; 
}

arma::mat element::getGlobalFrame(global_matrix matrix, int global_size) {
    if (matrix == global_matrix::k) {
        arma::mat global_frame(global_size, global_size, arma::fill::zeros);

        std::vector<int> gDims = { x1_g_, y1_g_, r1_g_, x2_g_, y2_g_, r2_g_ };
        for (int i = 0; i < gDims.size(); i++) {
            for (int j = 0; j < gDims.size(); j++) {
                global_frame(gDims[i] - 1, gDims[j] - 1) = global_k(i, j);
            }
        }
        return global_frame;
    }
    else {
        arma::mat global_frame(global_size, global_size, arma::fill::zeros);

        std::vector<int> gDims = { x1_g_, y1_g_, r1_g_, x2_g_, y2_g_, r2_g_ };
        for (int i = 0; i < gDims.size(); i++) {
            for (int j = 0; j < gDims.size(); j++) {
                global_frame(gDims[i] - 1, gDims[j] - 1) = global_m(i, j);
            }
        }
        return global_frame;
    }
}

// searches the frame structure and finds all coincident points 
// and numbers them with the same globals, xg, yg, rg
void frame::NumberSharedPoints(int xg, int yg, int rg, std::array<float, 2> shared_pt) {
    
    for (element* e : elements_) {
        for (int p = 0; p < 2; p++) {
            if (e->pts()[p] == shared_pt) 
                e->setGlobalxyr(xg, yg, rg, p);
        }
    }
}

void frame::CreateGCS() {
    // FindCoincendentPoints() (turn into function at some point )
    std::vector<std::array<float, 2>> copts = FindCoincidentPoints(); 

    // vector holding the coincident points that have been numbered 
    std::vector<std::array<float, 2>> copts_numbered;  

    // init globals
    int xg = 1, yg = 2, rg = 3; 

    // iterate through every element and every point in each element 
    for (element* e : elements_) {
        for (int p = 0; p < 2; p++) {
            // find out if element is a unique point 
            // will be unique if it is not found in coincident point vector, therefore: 
            if (std::find(copts.begin(), copts.end(), e->pts()[p]) == copts.end()) { // true = element is unique 
                e->setGlobalxyr(xg, yg, rg, p);
                xg += 3; 
                yg += 3; 
                rg += 3; 
            }
            else { // is coincident
                // test to see if it has been numbered already 
                // if current points is not found int copts_numbered vector, add it to it and increase the globals 
                if (std::find(copts_numbered.begin(), copts_numbered.end(), e->pts()[p]) == copts_numbered.end()) {
                    copts_numbered.push_back(e->pts()[p]);
                    e->setGlobalxyr(xg, yg, rg, p);
                    NumberSharedPoints(xg, yg, rg, e->pts()[p]);
                    xg += 3; 
                    yg += 3; 
                    rg += 3; 
                }
                else { /* do nothing, because it's already added and point has been numbered */ }
            }
        }
    }
    GCS_created = true; 
}

// iterates through all elements and finds the maximum GCS dim
int frame::GetMaximumGCSDim() {
    int max = 0; 
    for (element* e : elements_) 
        max = std::max(max, e->getMaxElementGlobalDim());
    return max;
}



void frame::global_k() {
#if DEBUG_GCS 
    printPointDataForStructure();
#endif
    CreateGCS();

#if DEBUG_GCS 
    printGCS();
#endif

    int global_size = GetMaximumGCSDim();
    arma::mat k_g(global_size, global_size, arma::fill::zeros); 
    arma::mat m_g(global_size, global_size, arma::fill::zeros); 
    

    // convert each globally transformed matrix of each element to the actual GLOBAL SIZED (12,12) we want
    for (element* e : elements_) {
        arma::mat e_k_g = e->getGlobalFrame(global_matrix::k, global_size);  // element global k matrix, sized for frame
        arma::mat e_m_g = e->getGlobalFrame(global_matrix::m, global_size);  // element global m matrix, sized for frame
        k_g += e_k_g; 
        m_g += e_m_g; 
    }
    k_g_ = k_g; 
    m_g_ = m_g; 

#if DEBUG_GCS
    printf("\nThe global stifness matrix is shown below:\n\n");
    std::cout << k_g_ << std::endl;
#endif

}