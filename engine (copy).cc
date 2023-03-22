#include "ini_configuration.h"
#include <fstream>
#include "l_parser.h"
#include "draw.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include "cmath"
#include "stack"
#include "Extra_functies.h"
using namespace std;

img::EasyImage generate_image(const ini::Configuration &configuration)
{
    img::EasyImage image;
    if(configuration["General"]["type"].as_string_or_die() == "2DLSystem"){
        vector<double> kleur = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
        vector<double> achtergrondkleur = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        int size = configuration["General"]["size"].as_int_or_die();
        LParser::LSystem2D l_system;
        ifstream fin(configuration["2DLSystem"]["inputfile"].as_string_or_die());
        fin>>l_system;
        fin.close();
        image= twoD(kleur,achtergrondkleur,l_system,size);
        return image;

    }
    string s=configuration["Figure0"]["type"].as_string_or_die();
    if(configuration["General"]["type"].as_string_or_die() == "Wireframe" xor configuration["General"]["type"].as_string_or_die() == "ZBufferedWireframe"  and configuration["Figure0"]["type"].as_string_or_die()=="LineDrawing"){
        vector<double> achtergrondkleur = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        vector<vector<double>>Points;
        vector<vector<double>>line;
        vector<vector<double>>center;
        double d=1;
        vector<vector<double>>kleur;
        vector<vector<vector<double>>>points;
        vector<vector<vector<double>>>lines;
        vector<double> rotateX;
        vector<double> rotateY;
        vector<double> rotateZ;
        vector<double> scale;
        double aantalF=configuration["General"]["nrFigures"].as_double_or_die();
        vector<double>eye=configuration["General"]["eye"].as_double_tuple_or_die();
        int size = configuration["General"]["size"].as_int_or_die();
        for (int i = 0; i <configuration["General"]["nrFigures"].as_double_or_die() ; ++i) {
            center.push_back(configuration["Figure" + to_string(i)]["center"].as_double_tuple_or_die());
            kleur.push_back(configuration["Figure" + to_string(i)]["color"].as_double_tuple_or_die());

            for (int j = 0; j <configuration["Figure" +to_string(i)]["nrPoints"].as_double_or_die() ; ++j) {
                Points.push_back(configuration["Figure" +to_string(i)]["point" + to_string(j)].as_double_tuple_or_die());
            }
            for (int k = 0; k <configuration["Figure" + to_string(i)]["nrLines"].as_double_or_die() ; ++k) {
                line.push_back(configuration["Figure" + to_string(i)]["line" + to_string(k)].as_double_tuple_or_die());
            }
             rotateX.push_back(configuration["Figure"+to_string(i)]["rotateX"].as_double_or_die()*(M_PI/180));
             rotateY.push_back(configuration["Figure"+to_string(i)]["rotateY"].as_double_or_die()*(M_PI/180));
             rotateZ.push_back(configuration["Figure"+to_string(i)]["rotateZ"].as_double_or_die()*(M_PI/180));
             scale.push_back(configuration["Figure"+to_string(i)]["scale"].as_double_or_die());
             points.push_back(Points);
             lines.push_back(line);
             Points.clear();
             line.clear();

        }

        image= threeD(kleur,achtergrondkleur,points,size,lines,center,eye,rotateX,rotateY,rotateZ,scale,aantalF,d);
        return image;


    }
    if(configuration["General"]["type"].as_string_or_die() == "Wireframe" xor configuration["General"]["type"].as_string_or_die() == "ZBufferedWireframe"  and configuration["Figure0"]["type"].as_string_or_die() != "3DLSystem"){
        vector<vector<double>>center;
        vector<string>type;
        vector<vector<double>>kleur;
        vector<double> rotateX;
        vector<double> rotateY;
        vector<double> rotateZ;
        vector<double> scale;
        double d=1;
        vector<int>n;
        vector<double>hoogte;
        vector<double>R;
        vector<double>r;
        vector<double>m;
        double aantalF=configuration["General"]["nrFigures"].as_double_or_die();
        vector<double> achtergrondkleur = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

        vector<double>eye=configuration["General"]["eye"].as_double_tuple_or_die();
        int size = configuration["General"]["size"].as_int_or_die();

        for (int i = 0; i <configuration["General"]["nrFigures"].as_double_or_die() ; ++i) {
            center.push_back(configuration["Figure" + to_string(i)]["center"].as_double_tuple_or_die());
            type.push_back(configuration["Figure" + to_string(i) ]["type"].as_string_or_die());
            kleur.push_back(configuration["Figure" + to_string(i)]["color"].as_double_tuple_or_die());
            rotateX.push_back(configuration["Figure"+to_string(i)]["rotateX"].as_double_or_die()*(M_PI/180));
            rotateY.push_back(configuration["Figure"+to_string(i)]["rotateY"].as_double_or_die()*(M_PI/180));
            rotateZ.push_back(configuration["Figure"+to_string(i)]["rotateZ"].as_double_or_die()*(M_PI/180));
            scale.push_back(configuration["Figure"+to_string(i)]["scale"].as_double_or_die());
            if (configuration["Figure"+to_string(i)]["n"].exists()){
                n.push_back(configuration["Figure"+to_string(i)]["n"].as_double_or_die());
            }
            if (!configuration["Figure"+to_string(i)]["n"].exists()){
                n.push_back(NULL);
            }
            if (configuration["Figure"+to_string(i)]["height"].exists()){
                hoogte.push_back(configuration["Figure"+to_string(i)]["height"].as_double_or_die());
            }
            if (!configuration["Figure"+to_string(i)]["height"].exists()){
                hoogte.push_back(NULL);
            }

            if (configuration["Figure"+to_string(i)]["m"].exists()){
                m.push_back(configuration["Figure"+to_string(i)]["m"].as_double_or_die());
            }
            if (!configuration["Figure"+to_string(i)]["m"].exists()){
                m.push_back(NULL);
            }
            if (configuration["Figure"+to_string(i)]["R"].exists()){
                R.push_back(configuration["Figure"+to_string(i)]["R"].as_double_or_die());
            }
            if (!configuration["Figure"+to_string(i)]["R"].exists()){
                R.push_back(NULL);
            }
            if (configuration["Figure"+to_string(i)]["r"].exists()){
                r.push_back(configuration["Figure"+to_string(i)]["r"].as_double_or_die());
            }
            if (!configuration["Figure"+to_string(i)]["r"].exists()){
                r.push_back(NULL);
            }
        }
         image = threeD1(kleur,center,rotateX,rotateY,rotateZ,scale,aantalF,eye,achtergrondkleur,type,size,n,hoogte,r,R,m,d);
        return image;

    }
    if(configuration["General"]["type"].as_string_or_die() == "Wireframe" xor configuration["General"]["type"].as_string_or_die() == "ZBufferedWireframe" and configuration["Figure0"]["type"].as_string_or_die()=="3DLSystem"){
        vector<vector<double>>center;
        vector<string>type;
        vector<vector<double>>kleur;
        vector<double> rotateX;
        vector<double> rotateY;
        vector<double> rotateZ;
        vector<double> scale;
        vector<int>n;
        vector<double>hoogte;
        vector<double>R;
        vector<double>r;
        vector<double>m;
        vector<double> achtergrondkleur = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        vector<double>eye=configuration["General"]["eye"].as_double_tuple_or_die();
        int size = configuration["General"]["size"].as_int_or_die();
        double aantalF=configuration["General"]["nrFigures"].as_double_or_die();
        for (int i = 0; i <aantalF ; ++i) {
            center.push_back(configuration["Figure" + to_string(i)]["center"].as_double_tuple_or_die());
            type.push_back(configuration["Figure" + to_string(i) ]["type"].as_string_or_die());
            kleur.push_back(configuration["Figure" + to_string(i)]["color"].as_double_tuple_or_die());
            rotateX.push_back(configuration["Figure"+to_string(i)]["rotateX"].as_double_or_die()*(M_PI/180));
            rotateY.push_back(configuration["Figure"+to_string(i)]["rotateY"].as_double_or_die()*(M_PI/180));
            rotateZ.push_back(configuration["Figure"+to_string(i)]["rotateZ"].as_double_or_die()*(M_PI/180));
            scale.push_back(configuration["Figure"+to_string(i)]["scale"].as_double_or_die());
            LParser::LSystem3D l_system;
            ifstream fin(configuration["Figure0"]["inputfile"].as_string_or_die());
            fin>>l_system;
            fin.close();
            image= threeDl(kleur,achtergrondkleur,l_system,size,rotateX,rotateY,rotateZ,center,eye,scale,aantalF);
    }
//        LParser::LSystem2D l_system;
//        ifstream fin(configuration["Figure0"]["inputfile"].as_string_or_die());
//        fin>>l_system;
//        fin.close();
//        image= threeDl(kleur,achtergrondkleur,l_system,size,rotateX,rotateY,rotateZ,center,eye,scale,aantalF);
        return image;
    }
    if(configuration["General"]["type"].as_string_or_die() == "Wireframe" xor configuration["General"]["type"].as_string_or_die() == "ZBuffering"  and configuration["Figure0"]["type"].as_string_or_die() != "3DLSystem"){
        vector<vector<double>>center;
        vector<string>type;
        vector<vector<double>>kleur;
        vector<double> rotateX;
        vector<double> rotateY;
        vector<double> rotateZ;
        vector<double> scale;
        double d=1;
        vector<int>n;
        vector<double>hoogte;
        vector<double>R;
        vector<double>r;
        vector<double>m;
        double aantalF=configuration["General"]["nrFigures"].as_double_or_die();
        vector<double> achtergrondkleur = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

        vector<double>eye=configuration["General"]["eye"].as_double_tuple_or_die();
        int size = configuration["General"]["size"].as_int_or_die();

        for (int i = 0; i <configuration["General"]["nrFigures"].as_double_or_die() ; ++i) {
            center.push_back(configuration["Figure" + to_string(i)]["center"].as_double_tuple_or_die());
            type.push_back(configuration["Figure" + to_string(i) ]["type"].as_string_or_die());
            kleur.push_back(configuration["Figure" + to_string(i)]["color"].as_double_tuple_or_die());
            rotateX.push_back(configuration["Figure"+to_string(i)]["rotateX"].as_double_or_die()*(M_PI/180));
            rotateY.push_back(configuration["Figure"+to_string(i)]["rotateY"].as_double_or_die()*(M_PI/180));
            rotateZ.push_back(configuration["Figure"+to_string(i)]["rotateZ"].as_double_or_die()*(M_PI/180));
            scale.push_back(configuration["Figure"+to_string(i)]["scale"].as_double_or_die());
            if (configuration["Figure"+to_string(i)]["n"].exists()){
                n.push_back(configuration["Figure"+to_string(i)]["n"].as_double_or_die());
            }
            if (!configuration["Figure"+to_string(i)]["n"].exists()){
                n.push_back(NULL);
            }
            if (configuration["Figure"+to_string(i)]["height"].exists()){
                hoogte.push_back(configuration["Figure"+to_string(i)]["height"].as_double_or_die());
            }
            if (!configuration["Figure"+to_string(i)]["height"].exists()){
                hoogte.push_back(NULL);
            }

            if (configuration["Figure"+to_string(i)]["m"].exists()){
                m.push_back(configuration["Figure"+to_string(i)]["m"].as_double_or_die());
            }
            if (!configuration["Figure"+to_string(i)]["m"].exists()){
                m.push_back(NULL);
            }
            if (configuration["Figure"+to_string(i)]["R"].exists()){
                R.push_back(configuration["Figure"+to_string(i)]["R"].as_double_or_die());
            }
            if (!configuration["Figure"+to_string(i)]["R"].exists()){
                R.push_back(NULL);
            }
            if (configuration["Figure"+to_string(i)]["r"].exists()){
                r.push_back(configuration["Figure"+to_string(i)]["r"].as_double_or_die());
            }
            if (!configuration["Figure"+to_string(i)]["r"].exists()){
                r.push_back(NULL);
            }
        }
        image = threeD2(kleur,center,rotateX,rotateY,rotateZ,scale,aantalF,eye,achtergrondkleur,type,size,n,hoogte,r,R,m,d);
        return image;

    }
    return image;
}

int main(int argc, char const* argv[])
{
        int retVal = 0;
        try
        {
                std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
                if (args.empty()) {
                        std::ifstream fileIn("filelist");
                        std::string filelistName;
                        while (std::getline(fileIn, filelistName)) {
                                args.push_back(filelistName);
                        }
                }
                for(std::string fileName : args)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(fileName);
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << fileName << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}
