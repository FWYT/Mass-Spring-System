#include "ParticleSystemLoader.h"
#include "GUILib/OBJReader.h"
#include <iostream>
#include <fstream>
#include "Utils/Logger.h"
#include <sstream>
#include <unordered_set>
#include "Constants.h"

struct SpringData {
	int p1;
	int p2;
	double stiffness;
};

struct ZeroLengthSpringData {
	int p;
	P3D x0;
	double stiffness;
};

ParticleSystem* ParticleSystemLoader::loadFromMSS(string filename) {
	ifstream infile(filename);
	string line;

	Logger::consolePrint("Opening %s", filename.c_str());

	if (!infile.is_open()) {
		Logger::consolePrint("Error opening file %s", filename.c_str());
		return 0;
	}

	vector<Particle> ps;
	vector<SpringData> sds;
	vector<ZeroLengthSpringData> zls;
	bool error = false;

	int lineNum = 0;

	while (getline(infile, line)) {
		lineNum++;
		istringstream iss(line);
		char c;
		if (line.at(0) == 'p') {
			double x1, x2, x3, v1, v2, v3, m;
			if (!(iss >> c >> x1 >> x2 >> x3 >> v1 >> v2 >> v3 >> m)) {
				Logger::consolePrint("Error reading particle on line %d", lineNum);
				error = true;
				break;
			}
			else {
				Logger::consolePrint("Particle with position (%f, %f, %f), velocity (%f, %f, %f), mass %f",
					x1, x2, x3, v1, v2, v3, m);
				Particle p;
				p.position = P3D(x1, x2, x3);
				p.velocity = V3D(v1, v2, v3);
				p.mass = m;
				ps.push_back(p);
			}
		}
		else if (line.at(0) == 's') {
			int p1, p2;
			double k;
			if (!(iss >> c >> p1 >> p2 >> k)) {
				Logger::consolePrint("Error reading spring on line %d", lineNum);
				error = true;
				break;
			}
			else {
				Logger::consolePrint("Spring between particles %d, %d, stiffness %f", p1, p2, k);
				SpringData sd;
				sd.p1 = p1;
				sd.p2 = p2;
				sd.stiffness = k;
				sds.push_back(sd);
			}
		}
		else if (line.at(0) == 'z') {
			int p;
			double x1, x2, x3, k;
			if (!(iss >> c >> p >> x1 >> x2 >> x3 >> k)) {
				Logger::consolePrint("Error reading zero-length spring on line %d", lineNum);
				error = true;
				break;
			}
			else {
				Logger::consolePrint("Zero-length spring pinning %d to (%f, %f, %f), stiffness %f",
					p, x1, x2, x3, k);
				ZeroLengthSpringData zl;
				zl.p = p;
				zl.x0 = P3D(x1, x2, x3);
				zl.stiffness = k;
				zls.push_back(zl);
			}
		}
		else {
			Logger::consolePrint("unrecognized");
		}
	}

	if (error) {
		Logger::consolePrint("Error reading file, no particles created.");
		return 0;
	}

	else {
		ParticleSystem* system = new ParticleSystem(ps);
		for (auto const &spring : sds) {
			system->addSpring(spring.p1, spring.p2, spring.stiffness);
		}
		for (auto const &spring : zls) {
			system->addZeroLengthSpring(spring.p, spring.x0, spring.stiffness);
		}
		return system;
	}
}

struct pair_hash {
	inline std::size_t operator()(const std::pair<int, int> & v) const {
		return v.first * 31 + v.second;
	}
};

ParticleSystem* ParticleSystemLoader::loadFromOBJ(string filename) {
	GLMesh* mesh = 0;
	try {
		mesh = OBJReader::loadOBJFile(filename.c_str());
	}

	catch (char* e) {
		Logger::consolePrint("%s", e);
		return 0;
	}
	
	vector<Particle> ps;
	unordered_set<pair<int, int>, pair_hash> edges;

	// First obtain vertices; these will become particles
	for (int i = 0; i < mesh->vertexCount; i++) {
		int base = 3 * i;
		Particle p;
		p.position = P3D(mesh->vertexList[base], mesh->vertexList[base + 1], mesh->vertexList[base + 2]);
		p.velocity = V3D(0, 0, 0);
		p.mass = DEFAULT_MASS;
		ps.push_back(p);
	}

	// Now obtain edges; these will become springs.
	// Need to iterate over triangles and add edges.
	for (auto &category : mesh->polygons->categories) {
		Logger::consolePrint("Sides per poly = %d", category->polyVertexCount);
		int numPolys = category->indexList.size() / category->polyVertexCount;
		Logger::consolePrint("Num polys = %d", numPolys);
		for (int i = 0; i < numPolys; i++) {
			int base = i * category->polyVertexCount;
			for (int j = 0; j < category->polyVertexCount; j++) {
				// Each pair of consecutive indices is another edge.
				int j_next = (j + 1) % category->polyVertexCount;
				pair<int, int> edge(category->indexList[base + j], category->indexList[base + j_next]);

				// If we haven't already added this edge, add it.
				if (edges.find(edge) == edges.end()) {
					edges.insert(edge);
				}
			}
		}
		/*
		Logger::consolePrint("Read triangle", edges.size());
		for (unsigned int i = 0; i < tri.indexes.size(); i++) {
			Logger::consolePrint("Read edge", edges.size());
			int next = (i + 1) % tri.indexes.size();

			pair<int, int> e(tri.indexes[i], tri.indexes[next]);
			
			// If we haven't already added this edge, add it.
			if (edges.find(e) == edges.end()) {
				edges.insert(e);
			}
		}*/
	}

	// Logger::consolePrint("%d edges", edges.size());

	ParticleSystem* system = new ParticleSystem(ps);
	for (auto &e : edges) {
		system->addSpring(e.first, e.second, DEFAULT_STIFFNESS);
	}
	system->setMesh(mesh);
	return system;
}