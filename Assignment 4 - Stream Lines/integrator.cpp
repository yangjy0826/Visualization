/*********************************************************************
*  Author  : Himangshu Saikia
*  Init    : Wednesday, September 20, 2017 - 12:04:15
*
*  Project : KTH Inviwo Modules
*
*  License : Follows the Inviwo BSD license model
*********************************************************************
*/

#include <labstreamlines/integrator.h>

namespace inviwo
{

	Integrator::Integrator()
	{
	}

	vec2 Integrator::sampleFromField(const VolumeRAM* vr, size3_t dims, const vec2& position)
	{
		// Sampled outside the domain!
		if (position[0] < 0 || position[0] > dims[0] - 1 || position[1] < 0 || position[1] > dims[1] - 1)
		{
			return vec2(0, 0);
		}

		int x0 = int(position[0]);
		int y0 = int(position[1]);

		// Leads to accessing only inside the volume
		// Coefficients computation takes care of using the correct values
		if (x0 == dims[0] - 1)
		{
			x0--;
		}
		if (y0 == dims[1] - 1)
		{
			y0--;
		}

		auto f00 = vr->getAsDVec2(size3_t(x0, y0, 0));
		auto f10 = vr->getAsDVec2(size3_t(x0 + 1, y0, 0));
		auto f01 = vr->getAsDVec2(size3_t(x0, y0 + 1, 0));
		auto f11 = vr->getAsDVec2(size3_t(x0 + 1, y0 + 1, 0));

		float x = position[0] - x0;
		float y = position[1] - y0;

		vec2 f;

		for (int i = 0; i < 2; i++)
		{
			f[i] = f00[i] * (1 - x) * (1 - y) + f01[i] * (1 - x) * y + f10[i] * x * (1 - y) + f11[i] * x * y;
		}

		return f;
	}


	// TODO: Implement a single integration step here

	vec2 Integrator::Euler(const VolumeRAM* vr, size3_t dims, const vec2& position, float stepSize)
	{
		//Access the vector field with sampleFromField(vr, dims, ...)
		auto v = sampleFromField(vr, dims, position);
		v *= stepSize;
		return position + v;
	}

	vec2 Integrator::Euler_direction(const VolumeRAM * vr, size3_t dims, const vec2 & position, float stepSize, int direction, bool inDirField)
	{
		auto v = sampleFromField(vr, dims, position);
		v *= direction;
		if (inDirField) {
			v *= 1. / sqrt(v[0] * v[0] + v[1] * v[1]);
		}
		v *= stepSize;
		return position + v;
	}

	//Stop the integration after a certain number of steps
	std::vector<vec2> Integrator::Euler_steps(const VolumeRAM* vr, size3_t dims, const vec2& position, float stepSize, int direction, bool inDirField, int steps)
	{
		std::vector<vec2> list;
		auto pt = position;
		list.push_back(pt);
		for (int i = 0; i < steps; i++) {
			if (direction == INTEGRATOR_DIRECTION_FORWARD) {
				pt = Euler_direction(vr, dims, pt, stepSize, +1, inDirField);
			}
			else if (direction == INTEGRATOR_DIRECTION_BACKWARD) {
				pt = Euler_direction(vr, dims, pt, stepSize, -1, inDirField);
			}
			list.push_back(pt);
		}
		return list;
	}

	//Stop the integration after a certain arc length of the stream line
	std::vector<vec2> Integrator::Euler_arc(const VolumeRAM* vr, size3_t dims, const vec2& position, float stepSize, int direction, bool inDirField, double arcLength)
	{
		std::vector<vec2> list;
		auto pt = position;
		auto pt_next = pt;
		list.push_back(pt);
		double totalPath = 0;
		do {
			if (direction == INTEGRATOR_DIRECTION_FORWARD) {
				pt_next = Euler_direction(vr, dims, pt, stepSize, +1, inDirField);
			}
			else if (direction == INTEGRATOR_DIRECTION_BACKWARD) {
				pt_next = Euler_direction(vr, dims, pt, stepSize, -1, inDirField);
			}
			auto diff = pt - pt_next;
			auto diff_mag = sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
			if (diff_mag < 1e-6) {
				break;
			}
			totalPath += diff_mag;
			list.push_back(pt_next);
			pt = pt_next;

		} while (totalPath < arcLength);
		return list;
	}

	//Stop the integration at the boundary of the domain
	std::vector<vec2> Integrator::Euler_boundary(const VolumeRAM* vr, size3_t dims, const vec2& position, float stepSize, int direction, bool inDirField, double minX, double  minY, double maxX, double maxY)
	{
		std::vector<vec2> list;
		auto pt = position;
		auto pt_next = pt;
		list.push_back(pt);
		double totalPath = 0;

		do {
			if (direction == INTEGRATOR_DIRECTION_FORWARD) {
				pt_next = Euler_direction(vr, dims, pt, stepSize, +1, inDirField);
			}
			else if (direction == INTEGRATOR_DIRECTION_BACKWARD) {
				pt_next = Euler_direction(vr, dims, pt, stepSize, -1, inDirField);
			}
			auto diff = pt - pt_next;
			auto diff_mag = sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
			if (diff_mag < 1e-6) {
				break;
			}
			if (pt.x<minX || pt.x>maxX || pt.y<minY || pt.y>maxY) {
				break;
			}
			totalPath += diff_mag;
			list.push_back(pt_next);
			pt = pt_next;

		} while (true);

		return list;
	}


	//Stop the integration at zeros of the vector field
	std::vector<vec2> Integrator::Euler_zeroes(const VolumeRAM* vr, size3_t dims, const vec2& position, float stepSize, int direction, bool inDirField)
	{
		std::vector<vec2> list;
		auto pt = position;
		auto pt_next = pt;
		list.push_back(pt);
		//double totalPath = 0;

		do {
			if (direction == INTEGRATOR_DIRECTION_FORWARD) {
				pt_next = Euler_direction(vr, dims, pt, stepSize, +1, inDirField);
			}
			else if (direction == INTEGRATOR_DIRECTION_BACKWARD) {
				pt_next = Euler_direction(vr, dims, pt, stepSize, -1, inDirField);
			}
			auto diff = pt - pt_next;
			auto diff_mag = sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
			if (diff_mag < 1e-6) {
				break;
			}
			//totalPath += diff_mag;
			list.push_back(pt_next);
			pt = pt_next;

		} while (true);

		return list;
	}

	//Stop the integration when the velocity becomes too slow
	std::vector<vec2> Integrator::Euler_slow(const VolumeRAM* vr, size3_t dims, const vec2& position, float stepSize, int direction, bool inDirField, double minSpeed)
	{
		std::vector<vec2> list;
		auto pt = position;
		auto pt_next = pt;
		list.push_back(pt);
		//double totalPath = 0;

		do {
			if (direction == INTEGRATOR_DIRECTION_FORWARD) {
				pt_next = Euler_direction(vr, dims, pt, stepSize, +1, inDirField);
			}
			else if (direction == INTEGRATOR_DIRECTION_BACKWARD) {
				pt_next = Euler_direction(vr, dims, pt, stepSize, -1, inDirField);
			}
			auto diff = pt - pt_next;
			diff *= 1. / stepSize;
			auto diff_mag = sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
			if (diff_mag < minSpeed + 1e-16) {
				break;
			}
			//totalPath += diff_mag;
			list.push_back(pt_next);
			pt = pt_next;

		} while (true);

		return list;
	}


	vec2 Integrator::RK4(const VolumeRAM* vr, size3_t dims, const vec2& position, float stepSize)
	{

		auto v1 = sampleFromField(vr, dims, position);
		auto x1 = position;
		auto v2 = sampleFromField(vr, dims, x1 + stepSize / 2 * v1);
		auto v3 = sampleFromField(vr, dims, x1 + stepSize / 2 * v2);
		auto v4 = sampleFromField(vr, dims, x1 + stepSize * v3);
		v1 *= 1. / 6;
		v2 *= 1. / 3;
		v3 *= 1. / 3;
		v4 *= 1. / 6;
		auto v = v1 + v2 + v3 + v4;
		v *= stepSize;
		return x1 + v;
	}

	// Allow integration in the direction field (normalized vector field)
	vec2 Integrator::RK4_direction(const VolumeRAM* vr, size3_t dims, const vec2& position, float stepSize, int direction, bool inDirField)
	{
		auto v1 = sampleFromField(vr, dims, position);
		auto x1 = position;
		auto v2 = sampleFromField(vr, dims, x1 + stepSize / 2 * v1);
		auto v3 = sampleFromField(vr, dims, x1 + stepSize / 2 * v2);
		auto v4 = sampleFromField(vr, dims, x1 + stepSize * v3);
		v1 *= 1. / 6;
		v2 *= 1. / 3;
		v3 *= 1. / 3;
		v4 *= 1. / 6;
		auto v = v1 + v2 + v3 + v4;
		v *= direction;
		if (inDirField) {
			auto M = sqrt(v[0] * v[0] + v[1] * v[1]);
			M = M > 1e-16 ? M : 1e-16;
			v *= 1. / M;
		}
		v *= stepSize;
		return x1 + v;
	}


	//Stop the integration after a certain number of steps
	std::vector<vec2> Integrator::RK4_steps(const VolumeRAM* vr, size3_t dims, const vec2& position, float stepSize, int direction, bool inDirField, int steps)
	{
		std::vector<vec2> list;
		auto pt = position;
		list.push_back(pt);
		for (int i = 0; i < steps; i++) {
			if (direction == INTEGRATOR_DIRECTION_FORWARD) {
				pt = RK4_direction(vr, dims, pt, stepSize, +1, inDirField);
			}
			else if (direction == INTEGRATOR_DIRECTION_BACKWARD) {
				pt = RK4_direction(vr, dims, pt, stepSize, -1, inDirField);
			}
			list.push_back(pt);
		}
		return list;
	}

	//Stop the integration after a certain arc length of the stream line
	std::vector<vec2> Integrator::RK4_arc(const VolumeRAM* vr, size3_t dims, const vec2& position, float stepSize, int direction, bool inDirField, double arcLength)
	{
		std::vector<vec2> list;
		auto pt = position;
		auto pt_next = pt;
		list.push_back(pt);
		double totalPath = 0;
		do {
			if (direction == INTEGRATOR_DIRECTION_FORWARD) {
				pt_next = RK4_direction(vr, dims, pt, stepSize, +1, inDirField);
			}
			else if (direction == INTEGRATOR_DIRECTION_BACKWARD) {
				pt_next = RK4_direction(vr, dims, pt, stepSize, -1, inDirField);
			}
			auto diff = pt - pt_next;
			auto diff_mag = sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
			if (diff_mag < 1e-6) {
				break;
			}
			totalPath += diff_mag;
			list.push_back(pt_next);
			pt = pt_next;

		} while (totalPath < arcLength);
		return list;
	}

	//Stop the integration at the boundary of the domain
	std::vector<vec2> Integrator::RK4_boundary(const VolumeRAM* vr, size3_t dims, const vec2& position, float stepSize, int direction, bool inDirField, double minX, double  minY, double maxX, double maxY)
	{
		std::vector<vec2> list;
		auto pt = position;
		auto pt_next = pt;
		list.push_back(pt);
		double totalPath = 0;

		do {
			if (direction == INTEGRATOR_DIRECTION_FORWARD) {
				pt_next = RK4_direction(vr, dims, pt, stepSize, +1, inDirField);
			}
			else if (direction == INTEGRATOR_DIRECTION_BACKWARD) {
				pt_next = RK4_direction(vr, dims, pt, stepSize, -1, inDirField);
			}
			auto diff = pt - pt_next;
			auto diff_mag = sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
			if (diff_mag < 1e-6) {
				break;
			}
			if (pt.x<minX || pt.x>maxX || pt.y<minY || pt.y>maxY) {
				break;
			}
			totalPath += diff_mag;
			list.push_back(pt_next);
			pt = pt_next;

		} while (true);

		return list;
	}


	//Stop the integration at zeros of the vector field
	std::vector<vec2> Integrator::RK4_zeroes(const VolumeRAM* vr, size3_t dims, const vec2& position, float stepSize, int direction, bool inDirField)
	{
		std::vector<vec2> list;
		auto pt = position;
		auto pt_next = pt;
		list.push_back(pt);
		//double totalPath = 0;

		do {
			if (direction == INTEGRATOR_DIRECTION_FORWARD) {
				pt_next = RK4_direction(vr, dims, pt, stepSize, +1, inDirField);
			}
			else if (direction == INTEGRATOR_DIRECTION_BACKWARD) {
				pt_next = RK4_direction(vr, dims, pt, stepSize, -1, inDirField);
			}
			auto diff = pt - pt_next;
			auto diff_mag = sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
			if (diff_mag < 1e-6) {
				break;
			}
			//totalPath += diff_mag;
			list.push_back(pt_next);
			pt = pt_next;

		} while (true);

		return list;
	}

	//Stop the integration when the velocity becomes too slow
	std::vector<vec2> Integrator::RK4_slow(const VolumeRAM* vr, size3_t dims, const vec2& position, float stepSize, int direction, bool inDirField, double minSpeed)
	{
		std::vector<vec2> list;
		auto pt = position;
		auto pt_next = pt;
		list.push_back(pt);
		//double totalPath = 0;

		do {
			if (direction == INTEGRATOR_DIRECTION_FORWARD) {
				pt_next = RK4_direction(vr, dims, pt, stepSize, +1, inDirField);
			}
			else if (direction == INTEGRATOR_DIRECTION_BACKWARD) {
				pt_next = RK4_direction(vr, dims, pt, stepSize, -1, inDirField);
			}
			auto diff = pt - pt_next;
			diff *= 1. / stepSize;
			auto diff_mag = sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
			if (diff_mag < minSpeed + 1e-16) {
				break;
			}
			//totalPath += diff_mag;
			list.push_back(pt_next);
			pt = pt_next;

		} while (true);

		return list;
	}

} // namespace

