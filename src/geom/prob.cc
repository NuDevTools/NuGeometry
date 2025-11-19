#include <iostream>
#include <typeinfo>
#include <climits>
#include <iomanip>
#include <typeinfo>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "geom/Vector3D.hh"
#include <numeric>
#include <random>
#include <map>
#include <string>
#include "geom/Volume.hh"
#include "geom/Ray.hh"
#include "geom/LineSegment.hh"
#include "geom/Parser.hh"
#include "geom/Random.hh"
#include "spdlog/spdlog.h"

using NuGeom::LogicalVolume;
using NuGeom::PhysicalVolume;
using NuGeom::Vector3D;


// class DetectorSim {
//     public:
//         void Setup(const std::string &geometry) {
//             // Create world from geometry

//             NuGeom::GDMLParser parser(geometry);
//             world = parser.GetWorld();
//         }

//         void SetupMaterials() {
//             int element_input_type;
//             std::cout<< "Do you want to 1) input the elements manually or 2) Import from a text file?" << std::endl;
//             std::cout<< "Choose 1 or 2" << std::endl;
//             std::cin>> element_input_type;

//             while(!std::cin.good() || ((element_input_type!= 1) && (element_input_type!= 2))){
//             std::cout<< "Invalid Input" << std::endl;
//             std::cin.clear();
//             std::cin.ignore(INT_MAX, '\n');

//             std::cout<< "1 or 2?";
//             std::cin>> std::setw(1) >> element_input_type;}


//             if (element_input_type==1){
//                 int num_of_elements;
//                 std::cout<< "How many different elements?" << std::endl; 
                
//             }

//         }

//         void SetupShapes() {
//             while(true) {
//                 std::string shape;
//                 std::cout << "What shape? (box, sphere, or exit) ";
//                 std::cin >> shape;
//                 if(shape == "box") {
//                     double x, y, z;
//                     std::cout << "Dimensions of box (x, y, z) ";
//                     std::cin >> x >> y >> z;
//                     // Check x, y, z > 0
//                     auto box = std::make_shared<NuGeom::Box>(x, y, z);
//                     shapes.push_back(box);
//                 } else if(shape == "sphere") {
//                     double r;
//                     std::cout << "Radius of sphere ";
//                     std::cin >> r;
//                     // Check r > 0
//                     auto sphere = std::make_shared<NuGeom::Sphere>(r);
//                     shapes.push_back(sphere);
//                 } else if(shape == "exit") {
//                     break;
//                 } else {
//                     std::cout << "Invalid shape!" << std::endl;
//                 }
//             }
//         }

//         void SetupVolumes() {
//             // 1. Select world shape
//             // 2. Ask if there is a shape inside
//             int depth = 1;
//             while(depth > 0) {
//             // 3. Select a shape 
//             // 4. Select a material
//             // 5. Rotate / Translate shape 
//             // 6. Add to world
//             // 7. Return to 2. until you say there is no more in world 
//             //    if yes: (depth++)
//             //    if no: (depth--)
//             }
//         }

//         std::vector<NuGeom::Material> GetMaterials() {
//             // return materials in detector
//         }

//         void GetMeanFreePath(const std::vector<double> &cross_section) {
//             // Fill result mfp
//             if(cross_section.size() != m_mats.size())
//                 throw "ERROR";

//             for(size_t i = 0; i < m_mats.size(); ++i) {
//                 m_mfp[m_mats[i]] = cross_section[i];
//             }
//         }

//         std::pair<Vector3D, NuGeom::Material> GetInteraction(const NuGeom::Ray &ray) {
//             auto segments = world.GetLineSegments(ray);
//             Vector3D point;
//             NuGeom::Material mat;
//             // Choose interaction point
//             return {point, mat};
//         }

//     private:
//         NuGeom::World world;
//         std::vector<std::shared_ptr<NuGeom::Shape>> shapes;
//         std::vector<NuGeom::Material> m_mats;
//         std::map<NuGeom::Material, double> m_mfp;

// };

NuGeom::Ray ShootRay(const NuGeom::Vector3D &corner1, const NuGeom::Vector3D &corner2) {
    auto [xmin, xmax] = std::minmax(corner1.X(), corner2.X());
    auto [ymin, ymax] = std::minmax(corner1.Y(), corner2.Y());
    auto [zmin, zmax] = std::minmax(corner1.Z(), corner2.Z());

    auto rand = NuGeom::Random::Instance();

    NuGeom::Vector3D position{rand.Uniform(xmin, xmax), rand.Uniform(ymin, ymax), rand.Uniform(zmin, zmax)};
    double costheta = rand.Uniform(0.0, 1.0);
    double sintheta = std::sqrt(1-costheta*costheta);
    double phi = rand.Uniform(0.0, 2*M_PI);
    NuGeom::Vector3D direction{sintheta*cos(phi), sintheta*sin(phi), costheta};

    return NuGeom::Ray(position, direction);
}

int main(){
    // Define materials in the detector
    NuGeom::Material mat("Water", 1.0, 2);
    NuGeom::Material mat1("Argon", 9.0, 1);
    mat.AddElement(NuGeom::Element("Hydrogen", 1, 1), 2);
    mat.AddElement(NuGeom::Element("Oxygen", 8, 16), 1);
    mat1.AddElement(NuGeom::Element("Argon", 18, 40), 1);
   
    // Define the interaction geometry
    // Define the inner detector
    auto inner_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{1, 1, 1}); // Define a 1x1x1 box
   
    auto inner_vol = std::make_shared<LogicalVolume>(mat, inner_box); 
    NuGeom::RotationX3D rot(45*M_PI/180.0);
    auto inner_pvol = std::make_shared<PhysicalVolume>(inner_vol, NuGeom::Transform3D{}, rot);


    // Define the outer detector
    auto outer_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{2, 2, 2});
    auto outer_vol = std::make_shared<LogicalVolume>(mat1, outer_box); 
    outer_vol->AddDaughter(inner_pvol);
    inner_vol->SetMother(outer_vol);
    NuGeom::RotationX3D rot2(30*M_PI/180.0);
    auto outer_pvol = std::make_shared<PhysicalVolume>(outer_vol, NuGeom::Transform3D{}, rot2);
    inner_pvol->SetMother(outer_pvol);

    // Define the "World"
    auto world_box = std::make_shared<NuGeom::Box>(NuGeom::Vector3D{4, 4, 4});
    auto world = std::make_shared<LogicalVolume>(mat, world_box);
    outer_vol->SetMother(world); 
    world -> AddDaughter(outer_pvol);

    // Calculate interaction location
    std::map<std::string,double> meanfreepaths{
        {"Water",1e+7},
        {"Air",1e+7},
        {"Argon",1e+7}
    };


    // Shoot Rays and get LineSegments
    double prob_max = 0;
    size_t ntrials = 1 << 22;
    for(size_t i = 0; i < ntrials; ++i) {
        NuGeom::Ray ray = ShootRay({-2, -2, -2}, {2, 2, -2});
        std::vector<NuGeom::LineSegment> segments0;
        world->GetLineSegments(ray, segments0);

        std::vector<double> probs0(segments0.size());
        std::vector<double> probs1(segments0.size());
        std::vector<double> seglength0(segments0.size());

        std::vector<std::string> material0(segments0.size());

        // Calculate probability to interact for each line segment
        for (size_t j=0; j<segments0.size(); ++j){
            material0[j]=segments0[j].GetMaterial().Name();

            // NOTE: This only works for l/meanfreepath tiny
            probs0[j]=segments0[j].Length()/meanfreepaths[material0[j]];
            // std::cout << -segments0[j].Length()/meanfreepaths[material0[j]] << std::endl;
            probs1[j]=1-exp(-segments0[j].Length()/meanfreepaths[material0[j]]);

            // For testing only:
            seglength0[j] = segments0[j].Length();
        }
        double normconst0=std::accumulate(probs0.begin(), probs0.end(), 0.0);
        if (normconst0>prob_max) {
            std::cout << "Trial: " << i << " Ray: O(" << ray.Origin() << ") D(" << ray.Direction() << ")\n";
            std::cout << "Hit " << segments0.size() << " segments\n";
            double normconst1=std::accumulate(probs1.begin(), probs1.end(), 0.0);
            std::cout << "Exact = " << normconst1 << " Approx = " << normconst0 << "\n";
            std::cout << "Old prob max = " << prob_max << "\n";
            std::cout << "New prob max = " << normconst0 << "\n";
            prob_max=normconst0;
        }
    }

    constexpr double safety_factor = 1.5;
    prob_max *= safety_factor;

    auto rand = NuGeom::Random::Instance();
    double nrays = 0;
    size_t nhits = 0;
    std::ofstream hist;
    hist.open("hit_locations.txt");

    // Send Ray for check to hit
    for(size_t i = 0; i < ntrials; ++i) {
        nrays += 1.0/prob_max;
        if(i % 10000 == 0) {
            std::cout << "Shot " << nrays << " rays\n";
        }
        NuGeom::Ray ray = ShootRay({-2, -2, -2}, {2, 2, -2});
        std::vector<NuGeom::LineSegment> segments0;

        world->GetLineSegments(ray, segments0);

        std::vector<double> probs0(segments0.size());
        std::vector<double> probs1(segments0.size());
        std::vector<double> seglength0(segments0.size());
        
        std::vector<std::string> material0(segments0.size());

        // Calculate probability to interact for each line segment
        for (size_t j=0; j<segments0.size(); ++j){
            material0[j]=segments0[j].GetMaterial().Name();
            // NOTE: This only works for l/meanfreepath tiny
            probs0[j]=segments0[j].Length()/meanfreepaths[material0[j]];
            // std::cout << -segments0[i].Length()/meanfreepaths[material0[i]] << std::endl;
            probs1[j]=1-exp(-segments0[j].Length()/meanfreepaths[material0[j]]);

            // For testing only:
            seglength0[j] = segments0[j].Length();
        }
        double normconst0=std::accumulate(probs0.begin(), probs0.end(), 0.0);
        // TODO: Determine if interaction occurs or not
        // 1. Define a warm-up step to find the maximum allowed total prob over a line segment (prob_max)
        // 2. Add some safety factor to scale by (prob_max = safety_factor * prob_max)
        // 3. Interaction occurs if prob / prob_max < R1 (R1 random number in [0, 1) )
        // 4. Determine where it occurs by throwing R2 and finding where \sum_i probs_i > R2 and pick region i
        // 5. Since probs << 1, select by choosing the point given by segement.start + (segement.end - segement.start) * (probs[i] - (R2 - sum_{j=0}^{i-1} probs[j]))/ probs[i]
        // 6. Return the interaction point to the user
        //
        // Setup Random number generator, seed, and generate random number in [0, 1)
        // auto rand = NuGeom::Random::Instance();
        // rand.Seed(2309182305);
        // double r1 = rand.Uniform(0.0, 1.0);
        if(normconst0 / prob_max < rand.Uniform(0.0, 1.0)) {
            continue;
        }
        nhits++;
        std::cout << nhits << " / " << nrays << "\n";

        for(size_t j=0; j<segments0.size(); ++j) probs0[j]/=normconst0;
        double sum = 0;
        size_t idx = 0;
        double r2 = rand.Uniform(0.0, 1.0);
        NuGeom::Vector3D position;
        for(size_t j = 0; j < probs0.size(); ++j) {
            sum += probs0[j]; 
            if(sum > r2) {
                idx = j;
                position = segments0[j].Start() + (segments0[j].End() - segments0[j].Start()) * (sum - r2)/probs0[j]; 
                break;
            }
        }

        std::cout << "Start: " << segments0[idx].Start() << "\n";
        std::cout << "End: " << segments0[idx].End() << "\n";
        std::cout << "Sum: " << sum << " R2: " << r2 << " Prob: " << probs0[idx] << " Dist: " << (sum - r2)/probs0[idx] << "\n";
        std::cout << "Hit Location: " << position << "\n";

        hist << position.X() << "," << position.Y() << "," << position.Z() << "\n";
        if(nhits == 10000) break;
    } 

    hist.close();

    return 0;

    /*

    std::random_device rd0;
    std::mt19937 rand_gen0(rd0());
    std::uniform_real_distribution<>dis0(0.0,1.0);

    size_t N0=10000;
    std::ofstream hist_sim;
    std::ofstream interpt0; 
    hist_sim.open("histogram_sim.txt");
    interpt0.open("interpt0.txt");

    size_t index0=0;
    Vector3D interaction_point0{};

    for (size_t i=0; i<N0; ++i){
        double rand0=dis0(rand_gen0);
        double sum_prob0=0;
        for (size_t j=0; j<probs0.size(); ++j){
            sum_prob0+=probs0[j];
            if (rand0<sum_prob0){
                index0=j;
                interaction_point0=segments0[j].Start()+(segments0[j].End()-segments0[j].Start())*((rand0+probs0[j]-sum_prob0)/(probs0[j]));
                break;}
        }
    interpt0<< interaction_point0.X() << "," << interaction_point0.Y() << "," << interaction_point0.Z() << "\n";
    hist_sim<< index0 <<"\n" ;
    } 
    hist_sim.close();
    interpt0.close();


// Sending out rays in random direction
    double N=200;
    std::ofstream interpt;
    std::ofstream interptXY;
    std::ofstream interptYZ;
    std::ofstream interptZX;
    interpt.open("interactionpt.txt");
    interptXY.open("XYinteraction.txt");
    interptYZ.open("YZinteraction.txt");
    interptZX.open("ZXinteraction.txt");

    std::random_device rd;
    std::mt19937 rand_gen(rd());
    std::uniform_real_distribution<>dis(0.0,1.0);

    // std::random_device rd2;
    // std::mt19937 rand_gen2(rd2());
    // std::uniform_real_distribution<>dis2(-1.0,1.0);
    


    for (int n=0; n<N; n++){               
        double rand=dis(rand_gen);
        double rand2=dis(rand_gen);
        double sum_prob=0;
            NuGeom::Ray ray({0, 0, -2}, {0, 1-(2*rand2), 1});
            std::vector<NuGeom::LineSegment> segments;
            world->GetLineSegments(ray, segments);

            std::vector<double> probs(segments.size());
            std::vector<double> seglength(segments.size());
            std::vector<double> invmeanfreepath(segments.size());
            std::vector<std::string> material(segments.size());

            for (size_t i=0; i<segments.size(); ++i){
                seglength[i]=segments[i].Length();
                material[i]=segments[i].GetMaterial().Name();
                if (material[i]=="Water")
                    {invmeanfreepath[i]=1/meanfreepaths["Water"];}
                else if (material[i]=="Air")
                    {invmeanfreepath[i]=1/meanfreepaths["Air"];}
                else if (material[i]=="Argon")
                    {invmeanfreepath[i]=1/meanfreepaths["Argon"];}

                //std::cout<< seglength[i] << " & " << 1-(2*rand) << "\n";
                // std::cout<< segments[i].X() << ' ' << segments[i].Y() << ' ' << segments[i].Z() << "\n";

                probs[i]=1.0-exp(-seglength[i]*invmeanfreepath[i]);}

                // for (size_t i=0; i<probs.size(); i++)
                // {
                //     std::cout<< probs[i] << "\n";
                // }

                double normconst= accumulate(probs.begin(), probs.end(), 0.0);
                
                // std::cout<< normconst << "\n"
                
                if(segments.size() > 5) {
                    std::cout << std::setprecision(16) << rand2 << std::endl;
                }

                for (size_t i=0; i<segments.size(); ++i)
                    {probs[i]=probs[i]/normconst;
                    //  std::cout<< probs[i] << "\n";
                    }


        Vector3D interaction_point{};
        for(size_t i=0; i<probs.size(); i++)
        {
            sum_prob+=probs[i];
            if(rand<sum_prob){ 
            interaction_point=segments[i].Start()+(segments[i].End()-segments[i].Start())*((rand+probs[i]-sum_prob)/(probs[i]));
            break;}
        }
        interpt<< interaction_point.X() << "," << interaction_point.Y() << "," << interaction_point.Z() << "\n";
        interptXY<< interaction_point.X() << "," << interaction_point.Y() << "\n";
        interptYZ<< interaction_point.Y() << "," << interaction_point.Z() << "\n";
        interptZX<< interaction_point.Z() << "," << interaction_point.X() << "\n";
    }
    interpt.close();
    interptXY.close();
    interptYZ.close();
    interptZX.close();

   
    return 0;
    
    
    //Now we are going to build a function that takes in the volumes, materials, ray direction etc from the user

    // Input of a map from material to mean free path 
    int map_input_type;
    std::map<std::string, double> input_mat_meanfreepathlen;

    std::cout<< "(1) Manually input the materials and mean free path, or (2) Import the Map from a txt file." << std::endl;
    std::cout<< "1 or 2?" << std::endl;
    //std::cin>> std::setw(1) >> map_input_type;

    while(!(std::cin>> std::setw(1) >> map_input_type) || ((map_input_type!= 1) && (map_input_type!= 2))){
        std::cout<< "Invalid Input. Choose 1 or 2." << std::endl;
        std::cin.clear();
        std::cin.ignore(INT_MAX, '\n');
    }
    

    if (map_input_type==1)
    {
        std::cout<< "How many different materials are there?" << std::endl;
        int number_of_mats;
        while(!(std::cin >> std::setw(1) >> number_of_mats) || (number_of_mats<=0)){
            std::cout<< "Error! Input has to be a positive integer. Try again!" << std::endl;
            std::cin.clear();
            std::cin.ignore(INT_MAX, '\n');
        }

        std::cin >> std::setw(0);

        std::string materialname;
        double meanfreepathlen;
        for (int i=0; i<number_of_mats; ++i){
            std::cout<< "What is the name of the material?" << std::endl;
            std::cin>> materialname;
            std::cout<< "What is the mean free path of the material?" << std::endl;
            std::cin>> meanfreepathlen;
            input_mat_meanfreepathlen.insert({materialname, meanfreepathlen});
        }
    } else if (map_input_type==2)
    {
        std::cout<< "Make sure the list is in the form: Name of material, mean free path ; and each mapping entry is in a new line." << std::endl; 
        std::cout<< "Print your file directory" << std::endl;
        std::string map_directory;
        std::cin>> map_directory;

        while(!std::cin.good()){
        std::cout<< "Invalid Input" << std::endl;
        std::cin.clear();
        std::cin.ignore(INT_MAX, '\n');
        std::cout<< "Print the directory";
        std::cin>> map_directory;
        }

        std::ifstream map_from_txt;
        map_from_txt.open(map_directory);
        //need to fix the directory issue
        if (map_from_txt.is_open()){
        std::string material_name;
        std::string line;
        std::string numstring;
        double meanfreepathlen;

        while(getline(map_from_txt, line))
        {
            std::stringstream ss(line);
            getline(ss, material_name, ',');
            getline(ss, numstring);
            meanfreepathlen = stod(numstring);
            input_mat_meanfreepathlen.insert({material_name, meanfreepathlen});
        }
        }
        map_from_txt.close();
    }

    //Input of a list of ray directions
    int ray_input_type;
    std::cout<< "(1) Manually input the rays, or (2) Import the rays from a txt file." << std::endl;
    std::cout<< "1 or 2?" << std::endl;
    std::cin>> std::setw(1) >> ray_input_type;

    while(!std::cin.good() || ((ray_input_type!= 1) && (ray_input_type!= 2))){
        std::cout<< "Invalid Input" << std::endl;
        std::cin.clear();
        std::cin.ignore(INT_MAX, '\n');

        std::cout<< "1 or 2?";
        std::cin>> std::setw(1) >> ray_input_type;
        }

    std::vector<NuGeom::Ray> vec_of_rays;

    if (ray_input_type==1){
        std::cin >> std::setw(0);
        int more_ray=1;
        while (more_ray==1)
        {
            std::cout<< "What are the x,y,z coordinates of the origin? (Type each and press enter)" << std::endl;
            double xo, yo, zo, xd, yd, zd;
            std::cin >> xo >> yo >> zo;
            //Vector3D originvec;
            std::cout<< "What are the x,y,z directions of the ray?" << std::endl;
            std::cin >> xd >> yd >> zd ;
            //Vector3D directionvec;
            NuGeom::Ray input_ray({xo, yo, zo}, {xd, yd, zd});
            vec_of_rays.push_back(input_ray);

            std::cout<< "Do you want to input more rays?-- (1) Yes, (2) No" << std::endl;
            int need_more_rays;
            std::cin>> need_more_rays;
            if (need_more_rays==1){more_ray=1;}
            else{more_ray=0;}
        }
    }
    else if (ray_input_type==2){
        std::cout<< "Make sure the list is in the form: x_origin, y_origin, z_origin, x_direction, y_direction, z_direction" << std::endl;
        std::cout<< "Print your file directory" << std::endl;
        std::string ray_directory;
        std::cin>> ray_directory;

        while(!std::cin.good()){
        std::cout<< "Invalid Input" << std::endl;
        std::cin.clear();
        std::cin.ignore(INT_MAX, '\n');
        std::cout<< "Print the directory";
        std::cin>> ray_directory;
        }

        std::ifstream ray_from_txt;
        ray_from_txt.open(ray_directory);

        if (ray_from_txt.is_open()){
        //std::cout<<"Open"<< std::endl;
        std::string vecline;
        std::string stringxo, stringyo, stringzo, stringxd, stringyd, stringzd;
        double xo, yo, zo, xd, yd, zd;

        while(getline(ray_from_txt, vecline))
        {
            std::stringstream ssv(vecline);
            getline(ssv, stringxo, ',');
            getline(ssv, stringyo, ',');
            getline(ssv, stringzo, ',');
            getline(ssv, stringxd, ',');
            getline(ssv, stringyd, ',');
            getline(ssv, stringzd);
            xo=stod(stringxo);
            yo=stod(stringyo);
            zo=stod(stringzo);
            xd=stod(stringxd);
            yd=stod(stringyd);
            zd=stod(stringzd);

            NuGeom::Ray input_ray({xo, yo, zo}, {xd, yd, zd});
            vec_of_rays.push_back(input_ray);
        }
        }
        ray_from_txt.close();
    }

    //Importing the objects
    //std::cout<< "What are the shapes and volumes of the objects?" << std::endl;
    //Import the shapes of the objects

    //map of cross sections to mean
    //std::map<double, double> map_crossection_meanfreepath;
    //code

    for (size_t i=0; i<vec_of_rays.size(); i++)
    {
        std::vector<NuGeom::LineSegment> segs;
        world->GetLineSegments(vec_of_rays[i],segs);
        std::cout << segs.size() << std::endl;

        std::vector<double> probabilities(segs.size());
        std::vector<double> segslength(segs.size());
        std::vector<double> meanfreepathsvec(segs.size());
        std::vector<std::string> materialsvec(segs.size());

        for (size_t j=0; j<segs.size(); ++j){
            segslength[j]=segs[j].Length();
            //std::cout<<segslength[j] << std::endl;
            materialsvec[j]=segs[j].GetMaterial().Name();
            //std::cout<< materialsvec[j] << std::endl;
            meanfreepathsvec[j]=input_mat_meanfreepathlen[materialsvec[j]];
            //std::cout<< meanfreepathsvec[j]<< std::endl;
            probabilities[j]=1.0-exp(-segslength[j]/meanfreepathsvec[j]);}

        double normalization_constant= accumulate(probabilities.begin(), probabilities.end(), 0.0);
        for (size_t k=0; k<segs.size(); ++k)
            {probabilities[k]=probabilities[k]/normalization_constant;
            std::cout<< probabilities[k] << " "<< typeid(probabilities[k]).name() << "\n";
            }
        }


    

    return 0;
    */
}
