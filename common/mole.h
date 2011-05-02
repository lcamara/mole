/*
 * Mole - Mobile Organic Localisation Engine
 * Copyright 2010 Nokia Corporation.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef MOLE_H_
#define MOLE_H_

#include <csignal>
#include <QtCore>
#include <QCoreApplication>
#include <QNetworkRequest>
#include <QNetworkReply>

#include "qglobal.h"
#define QT_FATAL_ASSERT 1

#define MOLE_ORGANIZATION "Nokia"
#define MOLE_APPLICATION  "Mole"
#define MOLE_VERSION      "27"

// d3f1ftrfyzkx8f.cloudfront.net <=> mole-static.research.nokia.com

// 184.73.196.5 <=> mole.research.nokia.com

#define MOLE_SERVER_URL    "http://200.141.73.238:8080"
#define MOLE_MAPSERVER_URL "http://200.141.73.238"
#define MOLE_MAPSERVER_URL_CHOP 12

// CHOP => number of characters between slash and end + 1 (this should be done programmatically....)

extern QSettings *settings;
extern QDir *app_root_dir;
extern bool debug;
extern bool verbose;

void init_mole_app(int argc, char *argv[], QCoreApplication *app, QString app_name, QString log_file_name);

void check_settings (bool reset);

QString get_current_version ();

//void start_timer (QTimer *timer, int period);
//void set_timer (QTimer *timer, bool increase, int &period, int min, int max);

extern const int UI_WIDTH;
extern const int UI_HEIGHT;

//extern const int UI_MAIN_WIDTH;
//extern const int UI_MAIN_HEIGHT;

//extern QString setting_username_key;
//extern QString setting_whereami_frequency_key;

#define setting_username_key "username"
#define setting_whereami_frequency_key "whereami_frequency"
#define setting_server_url_key "server_url"
#define setting_version_key "version"
#define setting_mapserver_url_key "mapserver_url"
#define setting_user_agent_key "user_agent"
#define setting_cookie_key "cookie"

extern QString setting_mole_server_url_value;
extern QString setting_server_url_value;
extern QString setting_mapserver_url_value;
extern QString setting_user_agent_value;

//extern QString setting_version_value;


extern QNetworkAccessManager *network_access_manager;

#endif /* MOLE_H_ */
