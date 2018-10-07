/*********************************************************************
*  Author  : Himangshu Saikia, Wiebke Koepp, ...
*  Init    : Monday, September 11, 2017 - 12:58:42
*
*  Project : KTH Inviwo Modules
*
*  License : Follows the Inviwo BSD license model
*********************************************************************
*/

#include <labmarchingsquares/marchingsquares.h>
#include <inviwo/core/util/utilities.h>

namespace inviwo
{

	// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
	const ProcessorInfo MarchingSquares::processorInfo_
	{
		"org.inviwo.MarchingSquares",      // Class identifier
		"Marching Squares",                // Display name
		"KTH Lab",                          // Category
		CodeState::Experimental,           // Code state
		Tags::None,                        // Tags
	};

	const ProcessorInfo MarchingSquares::getProcessorInfo() const
	{
		return processorInfo_;
	}


	MarchingSquares::MarchingSquares()
		:Processor()
		, inData("volumeIn")
		, meshOut("meshOut")
		, propShowGrid("showGrid", "Show Grid")
		, propDeciderType("deciderType", "Decider Type")
		//, propShowBothDeciders("propShowBothDeciders", "Show Both Deciders")
		, propMultiple("multiple", "Iso Levels")
		, propIsoValue("isovalue", "Iso Value")
		, propGridColor("gridColor", "Grid Lines Color", vec4(0.0f, 0.0f, 0.0f, 1.0f),
			vec4(0.0f), vec4(1.0f), vec4(0.1f),
			InvalidationLevel::InvalidOutput, PropertySemantics::Color)
		, propIsoColor("isoColor", "Color", vec4(0.0f, 0.0f, 1.0f, 1.0f),
			vec4(0.0f), vec4(1.0f), vec4(0.1f),
			InvalidationLevel::InvalidOutput, PropertySemantics::Color)
		, propIsoColorMin("isoColorMin", "ColorMin", vec4(0.0f, 0.0f, 1.0f, 1.0f),
			vec4(0.0f), vec4(1.0f), vec4(0.1f),
			InvalidationLevel::InvalidOutput, PropertySemantics::Color)
		, propIsoColorMax("isoColorMax", "ColorMax", vec4(1.0f, 0.0f, 0.0f, 1.0f),
			vec4(0.0f), vec4(1.0f), vec4(0.1f),
			InvalidationLevel::InvalidOutput, PropertySemantics::Color)
		, propNumContours("numContours", "Number of Contours", 1, 1, 50, 1)
		, propIsoTransferFunc("isoTransferFunc", "Colors", &inData)
	{
		// Register ports
		addPort(inData);
		addPort(meshOut);

		// Register properties
		addProperty(propShowGrid);
		addProperty(propGridColor);

		addProperty(propDeciderType);
		propDeciderType.addOption("midpoint", "Mid Point", 0);
		propDeciderType.addOption("asymptotic", "Asymptotic", 1);

		addProperty(propMultiple);

		propMultiple.addOption("single", "Single", 0);
		addProperty(propIsoValue);
		addProperty(propIsoColor);
		addProperty(propIsoColorMin);
		addProperty(propIsoColorMax);

		propMultiple.addOption("multiple", "Multiple", 1);
		addProperty(propNumContours);
		addProperty(propIsoTransferFunc);

		// The default transfer function has just two blue points
		propIsoTransferFunc.get().clearPoints();
		propIsoTransferFunc.get().addPoint(vec2(0.0f, 1.0f), propIsoColorMin.get());
		propIsoTransferFunc.get().addPoint(vec2(1.0f, 1.0f), propIsoColorMax.get());
		//propIsoTransferFunc.get().addPoint(vec2(0.0f, 1.0f), vec4(0.0f, 0.0f, 1.0f, 1.0f));
		//propIsoTransferFunc.get().addPoint(vec2(1.0f, 1.0f), vec4(0.0f, 0.0f, 1.0f, 1.0f));
		propIsoTransferFunc.setCurrentStateAsDefault();

		util::hide(propGridColor, propNumContours, propIsoTransferFunc);

		// Show the grid color property only if grid is actually displayed
		propShowGrid.onChange([this]()
		{
			if (propShowGrid.get())
			{
				util::show(propGridColor);
			}
			else
			{
				util::hide(propGridColor);
			}
		});

		// Show options based on display of one or multiple iso contours
		propMultiple.onChange([this]()
		{
			if (propMultiple.get() == 0)
			{
				util::show(propIsoValue, propIsoColor);
				util::hide(propNumContours, propIsoTransferFunc);
			}
			else
			{
				//util::hide(propIsoValue);
				//util::show(propIsoColor, propNumContours);

				// TODO (Bonus): Comment out above if you are using the transfer function
				// and comment in below instead
				util::hide(propIsoValue, propIsoColor);
				util::show(propNumContours, propIsoTransferFunc);
			}
		});

	}

	void MarchingSquares::process()
	{


		propIsoTransferFunc.get().clearPoints();
		propIsoTransferFunc.get().addPoint(vec2(0.0f, 1.0f), propIsoColorMin.get());
		propIsoTransferFunc.get().addPoint(vec2(1.0f, 1.0f), propIsoColorMax.get());
		//propIsoTransferFunc.get().addPoint(vec2(0.0f, 1.0f), vec4(0.0f, 0.0f, 1.0f, 1.0f));
		//propIsoTransferFunc.get().addPoint(vec2(1.0f, 1.0f), vec4(0.0f, 0.0f, 1.0f, 1.0f));
		propIsoTransferFunc.setCurrentStateAsDefault();

		if (!inData.hasData()) {
			return;
		}

		// This results in a shared pointer to a volume
		auto vol = inData.getData();

		// Extract the minimum and maximum value from the input data
		const double minValue = vol->dataMap_.valueRange[0];
		const double maxValue = vol->dataMap_.valueRange[1];

		// Set the range for the isovalue to that minimum and maximum
		propIsoValue.setMinValue(minValue);
		propIsoValue.setMaxValue(maxValue);

		// You can print to the Inviwo console with Log-commands:
		LogProcessorInfo("This scalar field contains values between " << minValue << " and " << maxValue << ".");
		// You can also inform about errors and warnings:
		// LogProcessorWarn("I am warning about something"); // Will print warning message in yellow
		// LogProcessorError("I am letting you know about an error"); // Will print error message in red
		// (There is also LogNetwork...() and just Log...(), these display a different source,
		// LogProcessor...() for example displays the name of the processor in the workspace while
		// Log...() displays the identifier of the processor (thus with multiple processors of the
		// same kind you would not know which one the information is coming from

		// Retreive data in a form that we can access it
		const VolumeRAM* vr = vol->getRepresentation< VolumeRAM >();
		const size3_t dims = vol->getDimensions();

		// Initialize mesh and vertices
		auto mesh = std::make_shared<BasicMesh>();
		std::vector<BasicMesh::Vertex> vertices;

		// Values within the input data are accessed by the function below
		// It's input is the VolumeRAM from above, the dimensions of the volume
		// and the indeces i and j of the position to be accessed where
		// i is in [0, dims.x-1] and j is in [0, dims.y-1]
		float valueat00 = getInputValue(vr, dims, 0, 0);
		LogProcessorInfo("Value at (0,0) is: " << valueat00);
		// You can assume that dims.z = 1 and do not need to consider others cases

		// TODO (Bonus) Gaussian filter
		// Our input is const, but you need to compute smoothed data and write it somewhere
		// Create an editable volume like this:
		//Volume volSmoothed(vol->getDimensions(), vol->getDataFormat());
		//auto vrSmoothed = volSmoothed.getEditableRepresentation<VolumeRAM>();
		// Values can be set with
		//vrSmoothed->setFromDouble(vec3(i,j,0), value);
		// getting values works with an editable volume as well
		//getInputValue(vrSmoothed, dims, 0, 0);


		// Gaussian 

		// true: apply gaussian filter
		// false: not apply gaussian filter
		//int GAUSS_FLAG = true;
		int GAUSS_FLAG = false;


		Volume volSmoothed(vol->getDimensions(), vol->getDataFormat());
		auto vrSmoothed = volSmoothed.getEditableRepresentation<VolumeRAM>();
		int radius = 3;
		float sigma = 0.5f;

		const VolumeRAM* m_vr = vr;

		if (GAUSS_FLAG == true)
			m_vr = gaussian_filter(vr, vrSmoothed, dims, sigma, radius);
		else
			m_vr = vr;


		// Grid

		// Properties are accessed with propertyName.get() 
		if (propShowGrid.get())
		{
			// TODO: Add grid lines of the given color 
			auto indexBufferGrid = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
			// draw row lines
			for (size_t j = 0; j < dims.y; j++) {
				vec2 v1 = vec2(0.0, float(j) / (dims.y - 1));
				vec2 v2 = vec2(1.0, float(j) / (dims.y - 1));
				drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid, vertices);
			}
			// draw colum lines 
			for (size_t i = 0; i < dims.x; i++) {
				vec2 v1 = vec2(float(i) / (dims.x - 1), 0.0);
				vec2 v2 = vec2(float(i) / (dims.x - 1), 1.0);
				drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid, vertices);
			}
			// The function drawLineSegments creates two vertices at the specified positions, 
			// that are placed into the Vertex vector defining our mesh. 
			// An index buffer specifies which of those vertices should be grouped into to make up lines/trianges/quads.
			// Here two vertices make up a line segment.
			//auto indexBufferGrid = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

			// Draw a line segment from v1 to v2 with a color, the coordinates in the final 
			// image range from 0 to 1 for both x and y
			//vec2 v1 = vec2(0.1, 0.1);
			//vec2 v2 = vec2(0.7, 0.7);
			//drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid, vertices);
		}

		// Iso contours

		if (propMultiple.get() == 0)
		{
			float isoValue = propIsoValue.get();
			LogProcessorInfo("IsoValue is: " << isoValue);
			std::vector<vec3> v;
			auto indexBufferGrid = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
			for (size_t j = 0; j < dims.y - 1; j++) {
				for (size_t i = 0; i < dims.x - 1; i++) {
					float value00 = getInputValue(m_vr, dims, i, j);
					float value01 = getInputValue(m_vr, dims, i, j + 1);
					float value10 = getInputValue(m_vr, dims, i + 1, j);
					float value11 = getInputValue(m_vr, dims, i + 1, j + 1);
					float min = getMinMax(value00, value01, value10, value11)[0];
					float max = getMinMax(value00, value01, value10, value11)[1];
					if (min < isoValue && max >= isoValue) {
						v = getIsoVertices(i, j, value00, value01, value10, value11, isoValue, dims);
						if (v.size() == 2) {
							vec2 v1 = v[0];
							vec2 v2 = v[1];
							drawLineSegment(v1, v2, propIsoColor.get(), indexBufferGrid, vertices);
							//LogProcessorInfo("color of IsoValue is: " << propIsoColor.get());
						}
						else if (v.size() == 4) {
							vec2 v1 = v[0];
							vec2 v2 = v[1];
							vec2 v3 = v[2];
							vec2 v4 = v[3];
							drawLineSegment(v1, v2, propIsoColor.get(), indexBufferGrid, vertices);
							drawLineSegment(v3, v4, propIsoColor.get(), indexBufferGrid, vertices);
							//LogProcessorInfo("IsoValue is: " << propIsoColor.get());
						}
					}
				}
			}
			// Midpoint Decider
			//if (propDeciderType.get() == 0) {

			//}
			// TODO: Draw a single isoline at the specified isovalue (propIsoValue) 
			// and color it with the specified color (propIsoColor)
		}
		else
		{
			// TODO: Draw the given number (propNumContours) of isolines between 
			// the minimum and maximum value
			int n = propNumContours.get();
			float step = float(maxValue - minValue) / (n + 1);
			//float step1 = float(maxValue - minValue) / (n + 1);
			std::vector<float> isoValues;
			std::vector<vec4> isoColors;
			for (int i = 1; i <= n; i++) {
				float f = minValue + step*i;
				isoValues.push_back(f);
				//isoColors.push_back(propIsoTransferFunc.get().sample((f - minValue) / (maxValue - minValue)));
				//isoColors.push_back(propIsoColor.get()); 
				//std::vector<vec4> c = propIsoTransferFunc.get().sample((f - minValue) / (maxValue - minValue));
				float c = (f - minValue) / (maxValue - minValue);
				//float c1=c;
				//if (n == 3) {
				float c1 = (c - (1.0 / (float)(n + 1))) / (1 - (2.0 / (float)(n + 1)));
				//}
				isoColors.push_back(propIsoTransferFunc.get().sample(c1));
				LogProcessorInfo("IsoValue is: " << f);
				//LogProcessorInfo("the color of IsoValue is: " << c);
				LogProcessorInfo("the color1 of IsoValue is: " << c1);

			}


			//for (int i = 0; i < isoColors.size(); i++) {

			//}
			//isoColors[isoColors.end()]-isoColors[isoColors.begin()];
			for (int iso = 0; iso < n; iso++) {
				float isoValue = isoValues[iso];
				std::vector<vec3> v;
				auto indexBufferGrid = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
				for (size_t j = 0; j < dims.y - 1; j++) {
					for (size_t i = 0; i < dims.x - 1; i++) {
						float value00 = getInputValue(m_vr, dims, i, j);
						float value01 = getInputValue(m_vr, dims, i, j + 1);
						float value10 = getInputValue(m_vr, dims, i + 1, j);
						float value11 = getInputValue(m_vr, dims, i + 1, j + 1);
						float min = getMinMax(value00, value01, value10, value11)[0];
						float max = getMinMax(value00, value01, value10, value11)[1];
						if (min < isoValue && max >= isoValue) {
							v = getIsoVertices(i, j, value00, value01, value10, value11, isoValue, dims);
							if (v.size() == 2) {
								vec2 v1 = v[0];
								vec2 v2 = v[1];
								drawLineSegment(v1, v2, isoColors[iso], indexBufferGrid, vertices);
								//LogProcessorInfo("the color of IsoValue is: " << isoColors[iso]);

							}
							else if (v.size() == 4) {
								vec2 v1 = v[0];
								vec2 v2 = v[1];
								vec2 v3 = v[2];
								vec2 v4 = v[3];
								drawLineSegment(v1, v2, isoColors[iso], indexBufferGrid, vertices);
								drawLineSegment(v3, v4, isoColors[iso], indexBufferGrid, vertices);
								//LogProcessorInfo("the color of IsoValue is: " << isoColors[iso]);

							}
						}
					}
				}
			}

			// TODO (Bonus): Use the transfer function property to assign a color
			// The transfer function normalizes the input data and sampling colors
			// from the transfer function assumes normalized input, that means
			// vec4 color = propIsoTransferFunc.get().sample(0.0f);
			// is the color for the minimum value in the data
			// vec4 color = propIsoTransferFunc.get().sample(1.0f);
			// is the color for the maximum value in the data

		}

		// Note: It is possible to add multiple index buffers to the same mesh,
		// thus you could for example add one for the grid lines and one for
		// each isoline
		// Also, consider to write helper functions to avoid code duplication
		// e.g. for the computation of a single iso contour

		mesh->addVertices(vertices);
		meshOut.setData(mesh);
	}

	double MarchingSquares::getInputValue(const VolumeRAM* data, const size3_t dims,
		const size_t i, const size_t j) {
		// Check if the indices are withing the dimensions of the volume
		if (i < dims.x && j < dims.y) {
			return data->getAsDouble(size3_t(i, j, 0));
		}
		else {
			LogProcessorError(
				"Attempting to access data outside the boundaries of the volume, value is set to 0");
			return 0;
		}
	}

	void MarchingSquares::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
		IndexBufferRAM* indexBuffer,
		std::vector<BasicMesh::Vertex>& vertices) {
		// Add first vertex
		indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
		// A vertex has a position, a normal, a texture coordinate and a color
		// we do not use normal or texture coordinate, but still have to specify them
		vertices.push_back({ vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color });
		// Add second vertex
		indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
		vertices.push_back({ vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color });
	}

	std::vector<float> MarchingSquares::getMinMax(float value00, float value01, float value10, float value11) {
		std::vector<float> vMinMax = { value00, value01, value10, value11 };
		//float min = std::numeric_limits<float>::max(), max = std::numeric_limits<float>::min();
		float min = vMinMax[0];
		float max = vMinMax[1];
		for (size_t i = 0; i < vMinMax.size(); i++) {
			min = vMinMax[i] < min ? vMinMax[i] : min;
			max = vMinMax[i] > max ? vMinMax[i] : max;
		}
		std::vector<float> ret;
		ret.push_back(min);
		ret.push_back(max);
		return ret;
	}
	std::vector<vec3> MarchingSquares::getIsoVertices(int xIdx, int yIdx, float value00, float value01, float value10, float value11, float isoValue, size3_t dims) {
		std::vector<vec3> vIso;
		vec3 isoPosBottom(((float)xIdx + (isoValue - value00) / (value10 - value00)) / (dims.x - 1), (float)yIdx / (dims.y - 1), 0);
		vec3 isoPosUp(vec3(((float)xIdx + (isoValue - value01) / (value11 - value01)) / (dims.x - 1), ((float)yIdx + 1) / (dims.y - 1), 0));
		vec3 isoPosLeft(vec3((float)xIdx / (dims.x - 1), ((float)yIdx + (isoValue - value00) / (value01 - value00)) / (dims.y - 1), 0));
		vec3 isoPosRight(vec3(((float)xIdx + 1) / (dims.x - 1), ((float)yIdx + (isoValue - value10) / (value11 - value10)) / (dims.y - 1), 0));

		//midpoint the average of the four
		float midPoint = (value00 + value01 + value10 + value11) / 4;
		//use formula on f(x_a,y_a) to get asymPoint
		//decision rule if c > f(x_a,y_a) connect (a,b) and (c,d)
		//else connect (a,d) and (b,c).
		float asymPoint = (value00*value11 - value10*value01) / (value11 + value00 - value10 - value01);


		std::vector<float> vList = { value00, value01, value10, value11,  isoValue };
		sort(vList.begin(), vList.end());
		//size_t isoValueIndex;
		size_t isoValueIndex = find(vList.begin(), vList.end(), isoValue) - vList.begin();

		// -+++ case or +++- case
		if ((isoValueIndex == 1) | (isoValueIndex == 3)) {
			if ((value00 == vList[0] && isoValueIndex == 1) | (value00 == vList[4] && isoValueIndex == 3)) {
				// value00 is - or +
				vIso.push_back(isoPosBottom);
				vIso.push_back(isoPosLeft);
			}
			else if ((value01 == vList[0] && isoValueIndex == 1) | (value01 == vList[4] && isoValueIndex == 3)) {
				// value01 is - or +
				vIso.push_back(isoPosUp);
				vIso.push_back(isoPosLeft);
			}
			else if ((value10 == vList[0] && isoValueIndex == 1) | (value10 == vList[4] && isoValueIndex == 3)) {
				// value10 is -
				vIso.push_back(isoPosBottom);
				vIso.push_back(isoPosRight);
			}
			else if ((value11 == vList[0] && isoValueIndex == 1) | (value11 == vList[4] && isoValueIndex == 3)) {
				// value11 is -
				vIso.push_back(isoPosUp);
				vIso.push_back(isoPosRight);
			}
		}
		// ++-- case
		else if (isoValueIndex == 2) {
			if ((value00 - isoValue)*(value10 - isoValue) > 0) {
				vIso.push_back(isoPosLeft);
				vIso.push_back(isoPosRight);
			}
			else if ((value00 - isoValue)*(value01 - isoValue) > 0) {
				vIso.push_back(isoPosUp);
				vIso.push_back(isoPosBottom);
			}
			else if ((value00 - isoValue)*(value11 - isoValue) > 0) {
				// Ambiguities
				if (propDeciderType.get() == 0) {
					// Midpoint Decider
					if ((midPoint - isoValue)*(value00 - isoValue) <0) {
						vIso.push_back(isoPosBottom);
						vIso.push_back(isoPosLeft);
						vIso.push_back(isoPosUp);
						vIso.push_back(isoPosRight);
					}
					else {
						vIso.push_back(isoPosBottom);
						vIso.push_back(isoPosRight);
						vIso.push_back(isoPosLeft);
						vIso.push_back(isoPosUp);
					}
				}
				else if (propDeciderType.get() == 1) {
					// Asymtotic Decider
					if ((asymPoint - isoValue)*(value00 - isoValue) <0) {
						vIso.push_back(isoPosBottom);
						vIso.push_back(isoPosLeft);
						vIso.push_back(isoPosUp);
						vIso.push_back(isoPosRight);
					}
					else {
						vIso.push_back(isoPosBottom);
						vIso.push_back(isoPosRight);
						vIso.push_back(isoPosLeft);
						vIso.push_back(isoPosUp);
					}
				}
			}
		}



		return vIso;
	}

	double MarchingSquares::gaussian(int x, int y, float sigma) {
		return (1.0 / (2.0 * sigma * M_PI) * exp(-(pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2))));
	}

	VolumeRAM* MarchingSquares::gaussian_filter(const VolumeRAM* vr, VolumeRAM* vrSmoothed, const size3_t dims, float sigma, const int radius)
	{
		for (int x = 0; x < dims.x; x++) {
			for (int y = 0; y < dims.y; y++) {

				double total = 0.0;
				double totalFilter = 0.0;

				int startX = (x - radius < 0) ? 0 : x - radius;
				int endX = (x + radius > dims.x - 1) ? dims.x - 1 : x + radius;
				int startY = (y - radius < 0) ? 0 : y - radius;
				int endY = (y + radius > dims.y - 1) ? dims.y - 1 : y + radius;

				for (int i = startX; i <= endX; i++) {
					for (int j = startY; j <= endY; j++) {
						double filterVal = gaussian(i - x, j - y, sigma);
						totalFilter += filterVal;
						total += getInputValue(vr, dims, i, j) * filterVal;
					}
				}
				vrSmoothed->setFromDouble(vec3(x, y, 0), total / totalFilter * 1.0);

			}
		}
		return vrSmoothed;
	}

} // namespace
