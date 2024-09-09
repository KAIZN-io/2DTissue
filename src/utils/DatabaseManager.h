#pragma once
#include "Struct.h"
#include <iostream>
#include <sqlite3.h>
#include <string>
#include <vector>

class DatabaseManager
{
  public:
    DatabaseManager(const std::string& db_name)
    {
        if (sqlite3_open(db_name.c_str(), &db) != SQLITE_OK)
        {
            std::cerr << "Failed to open database: " << sqlite3_errmsg(db) << std::endl;
            throw std::runtime_error("Failed to open SQLite database");
        }

        std::string createTableSystem = "CREATE TABLE IF NOT EXISTS system_info (order_parameter REAL);";
        std::string createTableParticles
            = "CREATE TABLE IF NOT EXISTS particles ("
              "x_UV REAL, y_UV REAL, x_velocity_UV REAL, y_velocity_UV REAL, "
              "alignment_UV REAL, x_3D REAL, y_3D REAL, z_3D REAL, neighbor_count INTEGER);";

        execSQL(createTableSystem);
        execSQL(createTableParticles);
    }

    ~DatabaseManager() { sqlite3_close(db); }

    void insertSystemInfo(const System& system)
    {
        std::string sql
            = "INSERT INTO system_info (order_parameter) VALUES (" + std::to_string(system.order_parameter) + ");";
        execSQL(sql);
    }

    void insertParticlesInfo(const std::vector<Particle>& particles)
    {
        for (const auto& particle : particles)
        {
            std::string sql = "INSERT INTO particles (x_UV, y_UV, x_velocity_UV, y_velocity_UV, alignment_UV, "
                              "x_3D, y_3D, z_3D, neighbor_count) VALUES ("
                              + std::to_string(particle.x_UV) + ", " + std::to_string(particle.y_UV) + ", "
                              + std::to_string(particle.x_velocity_UV) + ", " + std::to_string(particle.y_velocity_UV)
                              + ", " + std::to_string(particle.alignment_UV) + ", " + std::to_string(particle.x_3D)
                              + ", " + std::to_string(particle.y_3D) + ", " + std::to_string(particle.z_3D) + ", "
                              + std::to_string(particle.neighbor_count) + ");";
            execSQL(sql);
        }
    }

  private:
    sqlite3* db;

    void execSQL(const std::string& sql)
    {
        char* errMsg = nullptr;
        if (sqlite3_exec(db, sql.c_str(), nullptr, nullptr, &errMsg) != SQLITE_OK)
        {
            std::cerr << "Failed to execute query: " << errMsg << std::endl;
            sqlite3_free(errMsg);
            throw std::runtime_error("Failed to execute SQLite query");
        }
    }
};
