#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/Viewer.h>
#include <sstream>
#include <tuple>
#include <vector>

#include "GA.h"
#include "Skeleton.h"

using namespace Eigen;
using namespace std;

int main(int argc, char **argv)
{
	const Skeleton *curveSkeleton = new Skeleton(argv[1]); // High resolution: argv[1]
	const Skeleton *IKSkeleton = new Skeleton(argv[2]); // Low resolution: argv[2]
	GA *ga = new GA(curveSkeleton, IKSkeleton);

	// Embed the IK-skeleton onto the curveSkeleton. 
	Skeleton optSkeleton = ga->run();

	// Embed the IK-skeleton onto the curveSkeleton manually
	/*vector<int> pairs(IKSkeleton->vertices.size());

	Skeleton *optSkeleton = new Skeleton(IKSkeleton, curveSkeleton, pairs);*/

	//optSkeleton->saveAsCG("../IK_mocap09_01.cg");
	optSkeleton.saveAsCG(argv[3]);

	// ---------------------------- Primitives required to plot the IK skeleton -------------------------------
	/*MatrixXd IKSkelVertices;

	IKSkeleton->getSkeletonVertices(IKSkelVertices);
	vector<tuple<int, int>> IKSkelEdges = IKSkeleton->getSkeletonEdges();*/
	// --------------------------------------------------------------------------------------------------------

	// --------------------------- Primitives required to plot the curve skeleton -----------------------------
	MatrixXd curveSkelVertices;

	curveSkeleton->getSkeletonVertices(curveSkelVertices);
	vector<tuple<int, int>> curveSkelEdges = curveSkeleton->getSkeletonEdges();
	// --------------------------------------------------------------------------------------------------------

	// --------------------------- Primitives required to plot the optimum skeleton ---------------------------
	MatrixXd optSkelVertices;

	optSkeleton.getSkeletonVertices(optSkelVertices);
	vector<tuple<int, int>> optSkelEdges = optSkeleton.getSkeletonEdges();
	// --------------------------------------------------------------------------------------------------------

	igl::opengl::glfw::Viewer viewer;

	/*for(int i = 0 ; i < IKSkelEdges.size() ; i++)
	{
		int u = get<0>(IKSkelEdges[i]);
		int v = get<1>(IKSkelEdges[i]);

		viewer.data().add_edges(IKSkelVertices.row(u) / 30.0, IKSkelVertices.row(v) / 30.0, RowVector3d(1, 0, 0));
	}

	for(int i = 0 ; i < IKSkelVertices.rows() ; i++)
	{
		stringstream label;

		label << i;

		viewer.data().add_label(IKSkelVertices.row(i) / 30.0, label.str());
	}*/

	for(int i = 0 ; i < curveSkelEdges.size() ; i++)
	{
		int u = get<0>(curveSkelEdges[i]);
		int v = get<1>(curveSkelEdges[i]);

		viewer.data().add_edges(curveSkelVertices.row(u) / 30.0, curveSkelVertices.row(v) / 30.0, RowVector3d(1, 0, 0));
	}

	for(int i = 0 ; i < curveSkelVertices.rows() ; i++)
	{
		stringstream label;

		label << i;

		viewer.data().add_label(curveSkelVertices.row(i) / 30.0, label.str());
	}

	for(int i = 0 ; i < optSkelEdges.size() ; i++)
	{
		int u = get<0>(optSkelEdges[i]);
		int v = get<1>(optSkelEdges[i]);

		viewer.data().add_edges(optSkelVertices.row(u) / 30.0, optSkelVertices.row(v) / 30.0, RowVector3d(0, 1, 0));
	}

	for(int i = 0 ; i < optSkelVertices.rows() ; i++)
	{
		stringstream label;

		label << i;

		viewer.data().add_label(optSkelVertices.row(i) / 30.0, label.str());
	}

	viewer.data().show_labels = true;

	igl::opengl::glfw::imgui::ImGuiMenu menu;
  	
  	menu.callback_draw_viewer_window = [](){};
  
  	viewer.plugins.push_back(&menu);

  	viewer.launch();

	return 0;
}
