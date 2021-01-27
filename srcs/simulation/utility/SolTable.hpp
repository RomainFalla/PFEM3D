#pragma once
#ifndef SOLTABLE_HPP_INCLUDED
#define SOLTABLE_HPP_INCLUDED

#include <string>

#define SOL_ALL_SAFETIES_ON 1
#include <sol/sol.hpp>

class SolTable
{
    public:
        SolTable() = default;

        SolTable(const std::string& tableName, const sol::state& state):
        m_tableName(tableName)
        {
            sol::table table = state[tableName];
            if(!table.valid())
                throw std::runtime_error("Could not find table" + tableName + "inside sol::state!");

            m_tableInternal = table;
        }

        SolTable(const std::string& tableName, const SolTable& solTable):
        m_tableName(tableName)
        {
            sol::table table = solTable.m_tableInternal[tableName];
            if(!table.valid())
                throw std::runtime_error("Could not find table" + tableName + "inside sol::state!");

            m_tableInternal = table;
        }

        SolTable(sol::table table):
        m_tableName("anonymous")
        {
            m_tableInternal = table;
        }

        ~SolTable() = default;

        template<typename T, typename... Args>
        T call(const std::string& functionName, Args... args) const noexcept
        {
            sol::protected_function_result res = m_tableInternal[functionName](m_tableInternal, args...);

            return res.get<T>();
        }

        template<typename... Args>
        void checkCall(const std::string& functionName, Args... args) const
        {
            auto object = m_tableInternal[functionName](m_tableInternal, args...);
            if(!object.valid())
            {
                sol::error err = object;
                throw std::runtime_error("Cannot find method " + functionName + " inside table " + m_tableName + ": " +
                      std::string(err.what()) + "!");
            }
        }

        template<typename... Args>
        bool checkCallNoThrow(const std::string& functionName, Args... args) const
        {
            auto object = m_tableInternal[functionName](m_tableInternal, args...);
            return object.valid();
        }

        template<typename T>
        T get(const std::string& propertyName) const noexcept
        {
            auto object = m_tableInternal[propertyName];
            return object.get<T>();
        }

        template<typename T>
        T checkAndGet(const std::string& propertyName) const
        {
            //std::cout << propertyName << std::endl;
            sol::object object = m_tableInternal[propertyName];
            //std::cout << "Merde" << std::endl;
            if(!object.valid())
            {
                sol::error err = object.as<sol::error>();
                throw std::runtime_error("Property " + propertyName + " was not found in table " + m_tableName + ": " +
                                         std::string(err.what()) + "!");
            }

            return object.as<T>();
        }

        void for_each(std::function<void(sol::object /*key*/, sol::object /*value*/)> f) const
        {
            return m_tableInternal.for_each(f);
        }

    private:
        sol::table m_tableInternal;
        std::string m_tableName;
};

#endif // SOLTABLE_HPP_INCLUDED
