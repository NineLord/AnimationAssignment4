
#include "igl/opengl/glfw/renderer.h"
#include "tutorial/sandBox/inputManager.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

// Assignment 4 //
bool read_Meshes(igl::opengl::glfw::Viewer* viewer, string file) {
	ifstream conFile;
	string line;
	string model_name;
	int times;
	int selectedIndex = 0;
	conFile.open(file, ios::in);
	if (conFile.is_open()) {
		while (getline(conFile, line)) {
			viewer->load_mesh_from_file(line);
			viewer->load_mesh_from_file(line);
			break;
		}
		conFile.close();		
		return true;
	}
	else {
		cout << "Error: No configuration. \n" << endl;
		return false;
	}
}

void adjustModels(igl::opengl::glfw::Viewer* viewer) {

	for (auto i = 0; i < 2; i++) {
		viewer->selected_data_index = i;
		viewer->data().MyRotateX(2.5);
		viewer->data().MyRotateY(0.35);
		float val;
		i == 0 ? val = 1.5 : val = -1.5;
		viewer->data().Translate(Vector3f(val, 0, 0));
		viewer->data().velocity = 0.02;
		i == 0 ? val = -1 : val = 1;
		viewer->data().direction = Vector3f(val, 0, 0);
	}
	// Adjusting the "camera"
	viewer->MyTranslate(Vector3f(0, 0, -2));
	// viewer->TranslateInSystem(viewer->MakeTrans(), Vector3f(0, 0, -2), true);
}

// Assignment 3 //

//bool read_Meshes(igl::opengl::glfw::Viewer* viewer, string file, int amountOfCy) {
//	ifstream conFile;
//	string line;
//	string modelName;
//	int times;
//	int selectedIndex = 0;
//	conFile.open(file, ios::in);
//	if (conFile.is_open()) {
//		while (getline(conFile, line)) {
//			if (line.empty()) break;
//			modelName = line;
//			modelName.erase(modelName.find_last_of('.'));
//			modelName.erase(0, modelName.find_last_of('\\') + 1);
//			bool amISphere = !strcmp(&modelName[0], "sphere");
//			!amISphere ? times = amountOfCy : times = 1;
//			for (int i = 0; i < times; i++) {
//				viewer->load_mesh_from_file(line);
//				viewer->data_list[selectedIndex++].model = modelName;
//			}
//		}
//		conFile.close();
//
//		igl::opengl::ViewerData* curr = nullptr;
//		igl::opengl::ViewerData* prev = nullptr;
//
//		for (int i = 0; i < viewer->data_list.size(); i++) {
//			if (strcmp(&(viewer->data_list[i].model)[0], "sphere")) {
//				curr = &viewer->data_list[i];
//				curr->son = prev;
//				if (prev) {
//					prev->father = curr;
//				}
//				prev = curr;
//			}
//		}
//
//		return true;
//	}
//	else {
//		cout << "Error: No configuration. \n" << endl;
//		return false;
//	}
//}

//void adjustModels(igl::opengl::glfw::Viewer* viewer, int times) {
//	bool first = true;
//	int i;
//	float counter = 0;
//	float lenOfCy;
//	Eigen::Vector3d m;
//	Eigen::Vector3d M;
//	Eigen::MatrixXd axisPoints(6, 3);
//	Eigen::MatrixXi axisLines(5, 2);
//	igl::opengl::ViewerData* curr = nullptr;
//
//	for (i = 0; i < viewer->data_list.size(); i++) {
//		curr = &viewer->data_list[i];
//		curr->show_overlay_depth = false;
//
//		if (!(strcmp(&curr->model[0], "sphere"))) {
//			curr->MyTranslate(Eigen::Vector3f(5, 0, 0));
//		}
//		else {
//			m = curr->V.colwise().minCoeff();
//			M = curr->V.colwise().maxCoeff();
//			// Calculating only once, for the first cylinder.
//			if (first) {
//				axisPoints <<
//					0, m(1), 0,							// Zero
//					(M(0) + m(0)) / 2 + 1, m(1), 0,		// X
//					-((M(0) + m(0)) / 2 + 1), m(1), 0,	// - X
//					0, M(1) + 0.5, 0,					// Y
//					0, m(1), (M(2) + m(2)) / 2 + 1,		// Z
//					0, m(1), -((M(2) + m(2)) / 2 + 1);	// - Z
//
//				axisLines <<
//					0, 1,
//					0, 2,
//					0, 3,
//					0, 4,
//					0, 5;
//
//				lenOfCy = M(1) - m(1);
//				viewer->lengthOfArm = lenOfCy* times;
//				first = false;
//			}
//
//			curr->bottom << 0, m(1), 0;
//			curr->top << 0, M(1), 0;
//			curr->topF << 0, M(1), 0, 1;
//			curr->bottomF << 0, m(1), 0, 1;
//			curr->point_size = 10;
//			curr->line_width = 3;
//
//			curr->add_points(axisPoints, Eigen::RowVector3d(1, 0, 0));
//
//			for (unsigned j = 0;j < axisLines.rows(); ++j) {
//				curr->add_edges
//				(
//					axisPoints.row(axisLines(j, 0)),
//					axisPoints.row(axisLines(j, 1)),
//					Eigen::RowVector3d(0, 0, 1)
//				);
//			}
//
//			if (curr->son != nullptr)
//				curr->MyTranslate(Eigen::Vector3f(0, lenOfCy, 0));
//
//			curr->SetCenterOfRotation(curr->bottom);
//		}
//
//	}
//	// Adjusting the "camera"
//	viewer->TranslateInSystem(viewer->MakeTrans(), Eigen::Vector3f(-2.5, -2, -7), true);
//
//}

int main(int argc, char* argv[])
{
	Display* disp = new Display(1000, 800, "Wellcome");
	Renderer renderer;
	igl::opengl::glfw::Viewer viewer;
	// Assignment 4 //

	if (!(read_Meshes(&viewer, "configuration.txt"))) return 1;
	adjustModels(&viewer);

	// Assignment 3 //

	//int cyNum = 4;
	//if (!(read_Meshes(&viewer, "configuration.txt", cyNum))) return 1;
	//adjustModels(&viewer, cyNum);

	Init(*disp);
	renderer.init(&viewer);
	disp->SetRenderer(&renderer);
	disp->launch_rendering(true);
	delete disp;
}

