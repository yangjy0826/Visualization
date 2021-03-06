/*********************************************************************
*  Author  : Himangshu Saikia, Wiebke Koepp
*  Init    : Tuesday, September 19, 2017 - 15:08:24
*
*  Project : KTH Inviwo Modules
*
*  License : Follows the Inviwo BSD license model
*********************************************************************
*/

#include <labstreamlines/eulerrk4comparison.h>
#include <labstreamlines/integrator.h>
#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/interaction/events/mouseevent.h>

namespace inviwo
{

	// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
	const ProcessorInfo EulerRK4Comparison::processorInfo_
	{
		"org.inviwo.EulerRK4Comparison",      // Class identifier
		"Euler RK4 Comparison",               // Display name
		"KTH Lab",                        // Category
		CodeState::Experimental,              // Code state
		Tags::None,                           // Tags
	};

	const ProcessorInfo EulerRK4Comparison::getProcessorInfo() const
	{
		return processorInfo_;
	}

	EulerRK4Comparison::EulerRK4Comparison()
		:Processor()
		, outMesh("meshOut")
		, inData("inData")
		, propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(0), vec2(1024), vec2(0.5))
		// TODO: Initialize additional properties
		// propertyName("propertyIdentifier", "Display Name of the Propery", 
		// default value (optional), minimum value (optional), maximum value (optional), increment (optional));
		// propertyIdentifier cannot have spaces
		, numberOfIntegration("numberOfIntegration", "Number Of Integration", 100, 0, std::numeric_limits<int>::max())
		, stepSize("stepSize", "Step Size", 0.1, 0.0, std::numeric_limits<float>::max())
		, mouseMoveStart("mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
			MouseButton::Left, MouseState::Press | MouseState::Move)
		, propEulerColor("eulerColor", "Euler Color", vec4(1.0f, 0.0f, 0.0f, 1.0f),
			vec4(0.0f), vec4(1.0f), vec4(0.1f),
			InvalidationLevel::InvalidOutput, PropertySemantics::Color)
		, propRK4Color("rk4Color", "RK4 Color", vec4(0.0f, 1.0f, 0.0f, 1.0f),
			vec4(0.0f), vec4(1.0f), vec4(0.1f),
			InvalidationLevel::InvalidOutput, PropertySemantics::Color)
	{
		// Register Ports
		addPort(outMesh);
		addPort(inData);

		// Register Properties
		addProperty(propStartPoint);
		addProperty(mouseMoveStart);

		// TODO: Register additional properties
		addProperty(numberOfIntegration);
		addProperty(stepSize);
		addProperty(propEulerColor);
		addProperty(propRK4Color);
	}

	double EulerRK4Comparison::getInputValue(const VolumeRAM* data, size3_t dims, size_t x, size_t y)
	{
		if (x < dims.x && y < dims.y)
		{
			return data->getAsDouble(size3_t(x, y, 0));
		}
		else
		{
			LogProcessorError("Trying to access data ouside of the boundary, value is set as 0");
			return 0;
		}
	}

	void EulerRK4Comparison::drawPoint(vec3& pt, vec4& xColor, std::vector<BasicMesh::Vertex>& vertices, int& pointCounter, inviwo::IndexBufferRAM* bufferLine, inviwo::IndexBufferRAM* bufferPoint)
	{

		float xMargin = 0.1, yMargin = 0.1;
		double xScale = (1 - 2 * xMargin), yScale = (1 - 2 * yMargin);
		vec3 origin(xMargin, yMargin, 0);
		vec3 axisScaler(xScale, yScale, 0);

		//add first point
		auto plt_pt = vec3(pt.x / (dims.x - 1), pt.y / (dims.y - 1), 0);
		pt = origin + plt_pt * axisScaler;
		LogProcessorInfo("@" + std::to_string(pointCounter) + " [ " + std::to_string(pt[0]) + " , " + std::to_string(pt[1]) + " ] ");

		vertices.push_back({ pt, vec3(0), vec3(0), xColor });
		bufferPoint->add(static_cast<std::uint32_t>(pointCounter));
		bufferLine->add(static_cast<std::uint32_t>(pointCounter));
		pointCounter++;


	}

	void EulerRK4Comparison::eventMoveStart(Event* event)
	{
		auto mouseEvent = static_cast<MouseEvent*>(event);
		vec2 mousePos = mouseEvent->posNormalized();
		// Denormalize to volume dimensions
		mousePos.x *= dims.x - 1;
		//mousePos.y *= dims.y - 1;
		// Update starting point
		propStartPoint.set(mousePos);
		event->markAsUsed();
	}

	void EulerRK4Comparison::process()
	{
		// Get input
		if (!inData.hasData())
		{
			return;
		}
		auto vol = inData.getData();

		// Retreive data in a form that we can access it
		const VolumeRAM* vr = vol->getRepresentation< VolumeRAM >();
		dims = vol->getDimensions();
		// The start point should be inside the volume (set maximum to the upper right corner)
		propStartPoint.setMaxValue(vec2(dims.x - 1, dims.y - 1));

		// Initialize mesh, vertices and index buffers for the two streamlines and the 
		auto mesh = std::make_shared<BasicMesh>();
		std::vector<BasicMesh::Vertex> vertices;
		auto indexBufferEuler = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
		auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
		auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

		// Draw start point
		vec2 startPoint = propStartPoint.get();
		//vertices.push_back({ vec3(startPoint.x / (dims.x - 1), startPoint.y / (dims.y - 1), 0),
			//vec3(0), vec3(0), vec4(0, 0, 0, 1) });
		indexBufferPoints->add(static_cast<std::uint32_t>(0));
		indexBufferEuler->add(static_cast<std::uint32_t>(0));
		indexBufferRK->add(static_cast<std::uint32_t>(0));

		// TODO: Implement the Euler and Runge-Kutta of 4th order integration schemes
		// and then integrate forward for a specified number of integration steps and a given stepsize 
		// (these should be additional properties of the processor)

		//// data print
		//LogProcessorInfo(" --- Data::begin ---");
		//for (int i = 0; i < dims.x; i++) {
		//	for (int j = 0; j < dims.y; j++) {
		//		auto v = getInputValue(vr, dims, i, j);
		//		LogProcessorInfo(" [ " + std::to_string(i) + " , " + std::to_string(j) + " ] = " + std::to_string(v));
		//	}
		//}
		//LogProcessorInfo(" --- Data::end ---");
		//LogProcessorInfo(" ");
		int pointCounter = 0;


		// Integrator::Euler(vr, dims, startPoint, ...);
		// Integrator::Rk4(vr, dims, startPoint, ...);
		//Integrator::Euler(vr, dims, startPoint, numberOfIntegration, stepSize);
		auto pt = startPoint;
		Integrator integrator;
		auto plt_pt = vec3(pt, 0);
		drawPoint(plt_pt, propEulerColor.get(), vertices, pointCounter, indexBufferEuler, indexBufferPoints);
		LogProcessorInfo("--- Euler::begin ---");
		for (int i = 0; i < numberOfIntegration.get(); i++) {
			pt = integrator.Euler(vr, dims, pt, stepSize.get());
			plt_pt = vec3(pt, 0);
			LogProcessorInfo("@" + std::to_string(i) + " [ " + std::to_string(plt_pt[0]) + " , " + std::to_string(plt_pt[1]) + " ] ");
			drawPoint(plt_pt, propEulerColor.get(), vertices, pointCounter, indexBufferEuler, indexBufferPoints);
		}
		LogProcessorInfo("--- Euler::edn ---");
		LogProcessorInfo(" ");


		//Integrator::RK4(vr, dims, startPoint, numberOfIntegration, stepSize);
		pt = startPoint;
		//	Integrator integrator;
		plt_pt = vec3(pt, 0);
		drawPoint(plt_pt, propRK4Color.get(), vertices, pointCounter, indexBufferRK, indexBufferPoints);
		LogProcessorInfo("--- RK4::begin ---");
		for (int i = 0; i < numberOfIntegration.get(); i++) {
			pt = integrator.RK4(vr, dims, pt, stepSize.get());
			plt_pt = vec3(pt, 0);
			LogProcessorInfo("@" + std::to_string(i) + " [ " + std::to_string(plt_pt[0]) + " , " + std::to_string(plt_pt[1]) + " ] ");
			drawPoint(plt_pt, propRK4Color.get(), vertices, pointCounter, indexBufferRK, indexBufferPoints);
		}
		LogProcessorInfo("--- RK4::edn ---");
		LogProcessorInfo(" ");
		mesh->addVertices(vertices);
		outMesh.setData(mesh);
	}

} // namespace

