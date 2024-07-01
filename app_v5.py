import streamlit as st
import folium
from geopy.distance import geodesic
from math import radians, degrees, sin, cos, atan2
from streamlit_folium import folium_static
import numpy as np
import pandas as pd
import json
import time
import gpxpy
import mysql.connector
from streamlit_js_eval import get_geolocation
from streamlit_geolocation import streamlit_geolocation

MYSQL_HOST = st.secrets["MYSQL_HOST"]
MYSQL_PORT = st.secrets["MYSQL_PORT"]
MYSQL_USR = st.secrets["MYSQL_USR"]
MYSQL_PWD = st.secrets["MYSQL_PWD"]
MYSQL_SCHEMA = st.secrets["MYSQL_SCHEMA"]


def haversine(coord1, coord2):
    R = 6371  # Radius of the earth in km
    lat1, lon1 = coord1
    lat2, lon2 = coord2

    dlat = radians(lat2 - lat1)
    dlon = radians(lon2 - lon1)
    a = (sin(dlat / 2) * sin(dlat / 2) +
         cos(radians(lat1)) * cos(radians(lat2)) *
         sin(dlon / 2) * sin(dlon / 2))
    c = 2 * atan2(np.sqrt(a), np.sqrt(1 - a))
    d = 1000 * R * c  # Distance in km

    return d


def insert_into_database(coordinates):
    # conn = 0
    try:
        conn = mysql.connector.connect(
            host=MYSQL_HOST,
            port=MYSQL_PORT,
            user=MYSQL_USR,
            password=MYSQL_PWD,
            database=MYSQL_SCHEMA
        )
        cursor = conn.cursor()
        query = "INSERT INTO `coursemarks` (`type`, `latitude`, `longitude`, `dropTimeUtc`) VALUES (%s, %s, %s, %s)"
        current_time = time.strftime('%Y-%m-%d %H:%M:%S')

        # Insert mark coordinates
        mark_types = ['RC', 'Pin', 'WG1', 'WG2']
        for i, mark in enumerate(coordinates['marks']):
            cursor.execute(
                query, (mark_types[i], mark[0], mark[1], current_time))

        # Insert boundary coordinates
        boundary_types = ['B1', 'B2', 'B3', 'B4']
        for i, boundary in enumerate(coordinates['boundary']):
            cursor.execute(
                query, (boundary_types[i], boundary[0], boundary[1], current_time))

        conn.commit()
        st.success('Coordinates uploaded to the database successfully!')
    except mysql.connector.Error as err:
        st.error(f"Error: {err}")
    finally:
        if conn.is_connected():
            cursor.close()
            conn.close()


def add_marker(map_obj, location, tooltip):
    folium.Marker(location=location, tooltip=tooltip).add_to(map_obj)


def move_point(origin, bearing, distance_nm):
    """Move a point by a given distance (in nautical miles) at a given bearing."""
    destination = geodesic(kilometers=distance_nm *
                           1.852).destination(origin, bearing)
    return [destination.latitude, destination.longitude]


def calculate_heading(coord1, coord2):
    """Calculate the heading between two coordinates."""
    lat1, lon1 = map(radians, coord1)
    lat2, lon2 = map(radians, coord2)
    dLon = lon2 - lon1
    y = sin(dLon) * cos(lat2)
    x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon)
    heading = (degrees(atan2(y, x)) + 360) % 360
    return heading


def main():
    st.title('Sailing Race Course Setup')

    # Slider for TWD
    twd = st.slider('True Wind Direction (TWD)', 0, 360, 0)

    location = streamlit_geolocation()

    st.write(location)

    if st.checkbox("Check my location"):
        loc = get_geolocation()
        st.write(f"Your coordinates are {loc}")

    # Dropdown menu for selecting input method
    input_method = st.sidebar.selectbox(
        'Select Input Method', ('Actual Position', 'Manual coordinates', 'GPX coordinates'))

    rc, pin, wg1, wg2 = [None, None], [None, None], [None, None], [None, None]
    waypoints = []

    if input_method == 'Manual coordinates':
        st.sidebar.header('Input Coordinates')

        # Input coordinates for the marks
        rc_lat = st.sidebar.text_input('RC Latitude', '41.4269294331907')
        rc_lon = st.sidebar.text_input('RC Longitude', '2.254188602965819')
        pin_lat = st.sidebar.text_input('Pin Latitude', '41.428403727639264')
        pin_lon = st.sidebar.text_input('Pin Longitude', '2.2526792367673165')
        wg1_lat = st.sidebar.text_input('WG1 Latitude', '41.43615599027518')
        wg1_lon = st.sidebar.text_input('WG1 Longitude', '2.2860348853217105')
        wg2_lat = st.sidebar.text_input('WG2 Latitude', '41.434486794864625')
        wg2_lon = st.sidebar.text_input('WG2 Longitude', '2.2865501736423925')

        rc = [float(rc_lat), float(rc_lon)]
        pin = [float(pin_lat), float(pin_lon)]
        wg1 = [float(wg1_lat), float(wg1_lon)]
        wg2 = [float(wg2_lat), float(wg2_lon)]
        
    elif input_method == 'Actual Position':
        st.sidebar.header('Ping the mark where you are')
        if st.sidebar.button("RC position"):
            #rc = [streamlit_geolocation(key="rc")['latitude'], streamlit_geolocation(key="rc")['longitude']]
            location = streamlit_geolocation()

            st.write(location)
            #loc = get_geolocation()
            #st.write(f"Your coordinates are {loc}")
            #rc = [streamlit_geolocation()['latitude'], streamlit_geolocation()['longitude']]
            #st.write(f"RC coordinates are {rc}")
            
        elif  st.sidebar.button("Pin Position"):
            def get_actual_position():
                
            pin = [streamlit_geolocation(key="pin")['latitude'], streamlit_geolocation(key="pin")['longitude']]
            st.write(f"Pin coordinates are {pin}")
            
        elif  st.sidebar.button("WG1 Position"):
            wg1 = [streamlit_geolocation(key="wg1")['latitude'], streamlit_geolocation(key="wg1")['longitude']]
            st.write(f"WG1 coordinates are {wg1}")
            
        elif  st.sidebar.button("WG2 Position"):
            wg2 = [streamlit_geolocation(key="wg2")['latitude'], streamlit_geolocation(key="wg2")['longitude']] #[get_geolocation()[0], get_geolocation()[1]]
            st.write(f"WG2 coordinates are {wg2}")
            
        else :
            st.sidebar.write("don't forget to ping all the marks mate!")
        

    elif input_method == 'GPX coordinates':
        st.sidebar.header('Upload GPX Files')
        rc_file = st.sidebar.file_uploader('Upload RC GPX file', type=['gpx'])
        pin_file = st.sidebar.file_uploader(
            'Upload Pin GPX file', type=['gpx'])
        wg1_file = st.sidebar.file_uploader(
            'Upload WG1 GPX file', type=['gpx'])
        wg2_file = st.sidebar.file_uploader(
            'Upload WG2 GPX file', type=['gpx'])

        def extract_waypoint(file):
            if file:
                gpx = gpxpy.parse(file)
                waypoint = gpx.waypoints[0]
                return [waypoint.latitude, waypoint.longitude]
            return [None, None]

        rc = extract_waypoint(rc_file)
        pin = extract_waypoint(pin_file)
        wg1 = extract_waypoint(wg1_file)
        wg2 = extract_waypoint(wg2_file)

    if None in rc or None in pin or None in wg1 or None in wg2:
        st.warning('Please provide coordinates for all marks.')
        return

    # Calculate race course axis
    course_heading = calculate_heading(rc, wg1)
    perpendicular_heading = (course_heading + 90) % 360
    perpendicular_heading_lee_gate = (calculate_heading(rc, pin) + 90) % 360
    perpendicular_heading_win_gate = (calculate_heading(wg2, wg1) + 90) % 360
    distance_start = haversine(rc, pin)
    distance_uwgate = haversine(wg1, wg2)
    start_bias = distance_start * \
        np.tan((perpendicular_heading_lee_gate-twd)*np.pi/180)
    winward_bias = distance_start * \
        np.tan((perpendicular_heading_win_gate-twd)*np.pi/180)


    # Slider for boundary width and length
    boundary_width = st.sidebar.slider('Boundary Width (NM)', 0.0, 2.0, 0.6)
    boundary_length = st.sidebar.slider('Boundary Length (NM)', 0.0, 2.0, 1.1)

    # Calculate boundary coordinates around the marks
    def calculate_boundaries(marks, boundary_width, boundary_length):
        rc, pin, wg1, wg2 = marks

        # Calculate the center points of the sides
        center_bottom = [(rc[0] + pin[0]) / 2, (rc[1] + pin[1]) / 2]
        center_top = [(wg1[0] + wg2[0]) / 2, (wg1[1] + wg2[1]) / 2]

        # Move these center points to create the boundary
        boundary_bottom_left = move_point(
            center_bottom, perpendicular_heading + 180, boundary_width / 2)
        boundary_bottom_right = move_point(
            center_bottom, perpendicular_heading, boundary_width / 2)
        boundary_top_left = move_point(
            center_top, perpendicular_heading + 180, boundary_width / 2)
        boundary_top_right = move_point(
            center_top, perpendicular_heading, boundary_width / 2)

        boundary_bottom_left = move_point(
            boundary_bottom_left, course_heading + 180, boundary_length / 2)
        boundary_bottom_right = move_point(
            boundary_bottom_right, course_heading + 180, boundary_length / 2)
        boundary_top_left = move_point(
            boundary_top_left, course_heading, boundary_length / 2)
        boundary_top_right = move_point(
            boundary_top_right, course_heading, boundary_length / 2)

        #
        return [boundary_bottom_left, boundary_bottom_right, boundary_top_right, boundary_top_left, boundary_bottom_left]

    marks = [rc, pin, wg1, wg2]
    boundary_coords = calculate_boundaries(
        marks, boundary_width, boundary_length)

    # Display map
    center = (np.array(rc) + np.array(wg1) + np.array(wg2) + np.array(pin)) / 4

    m = folium.Map(location=center, zoom_start=14)
    marks_data = pd.DataFrame(columns=['Latitude', 'Longitude', 'name'], index=[
                              'rc', 'pin', 'Wg1', 'Wg2'])
    marks_data.loc['rc'] = [rc[0], rc[1], 'RC']
    marks_data.loc['pin'] = [pin[0], pin[1], 'Pin']
    marks_data.loc['Wg1'] = [wg1[0], wg1[1], 'WG1']
    marks_data.loc['Wg2'] = [wg2[0], wg2[1], 'WG2']
    for i in range(0, len(marks_data)):
        folium.Marker(
            location=[marks_data.iloc[i]['Latitude'],
                      marks_data.iloc[i]['Longitude']],
            popup=marks_data.iloc[i]['name'],
            icon=folium.DivIcon(f"""
            <div style="color: blue">{marks_data.iloc[i]['name']}</div>
                """)
        ).add_to(m)

    add_marker(m, rc, 'RC')
    add_marker(m, pin, 'Pin')
    add_marker(m, wg1, 'WG1')
    add_marker(m, wg2, 'WG2')

    folium.PolyLine(locations=boundary_coords, color='red',
                    dash_array='5, 10').add_to(m)

    st.subheader('Race Course Details')
    st.write(f'Course Axis Heading: {course_heading:.2f}°')
    #st.write(f'Leeward gate square at: {perpendicular_heading_lee_gate:.2f}°')
    #st.write(f'Leeward gate distance: {distance_start:.2f}m')
    #st.write(f'Leeward gate bias: {start_bias:.2f}m')

    #st.write(f'Windward gate square at: {perpendicular_heading_win_gate:.2f}°')

    #st.write(f'Windward gate distance: {distance_uwgate:.2f}m')
    #st.write(f'Windward gate bias: {winward_bias:.2f}m')
    
    recap_table = pd.DataFrame(columns=['Leeward Gate','Winward Gate'], index=['Gate square at', 'Gate distance', 'Gate bias'])
    recap_table.loc['Gate square at', 'Leeward Gate'] = f'{perpendicular_heading_lee_gate:.2f}°'
    recap_table.loc['Gate distance', 'Leeward Gate'] = f'{distance_start:.2f}m'
    recap_table.loc['Gate bias', 'Leeward Gate'] = f'{start_bias:.2f}m'

    recap_table.loc['Gate square at', 'Winward Gate'] = f'{perpendicular_heading_win_gate:.2f}°'
    recap_table.loc['Gate distance', 'Winward Gate'] = f'{distance_uwgate:.2f}m'
    recap_table.loc['Gate bias', 'Winward Gate'] = f'{winward_bias:.2f}m'
    
    st.dataframe(recap_table)
    
    folium_static(m)

    # Create JSON structure
    json_data = {
        "Boundary": boundary_coords,
        "RaceArea": [],
        "OtherAreas": [],
        "CreationTimeDate": time.time(),
        "CreationTime": time.time(),
        "RaceStartTime": time.time() + 600,
        "EntryBeforeStartBlue": 130.0,
        "EntryBeforeStartYellow": 120.0,
        "RaceType": "Match",
        "Participants": ["RaceBoat", "Ghost"],
        "OtherAssets": [],
        "CourseSequence": [
            {
                "Name": "0",
                "Rounding": "",
                "Gate": [
                    {"Name": "StartBoat", "TargetLocation": rc,
                        "BoatID": "StartBoat", "ForceTarget": False},
                    {"Name": "Pin", "TargetLocation": pin,
                        "BoatID": "Pin", "ForceTarget": False}
                ]
            },
            {
                "Name": "1",
                "Rounding": "",
                "Gate": [
                    {"Name": "TopMarkPort", "TargetLocation": wg1,
                        "BoatID": "TopMarkPort", "ForceTarget": False},
                    {"Name": "TopMarkStbd", "TargetLocation": wg2,
                        "BoatID": "TopMarkStbd", "ForceTarget": False}
                ]
            }
        ],
        "CourseDimensions": {
            "CourseWidth": boundary_width,
            "StartBoxDepth": boundary_length,
            "CourseLength": boundary_length,
            "CourseAxis": course_heading,
            "MarkZoneSize": 60.0,
            "BoundaryZoneSize": 115.0
        }
    }

    # Download JSON button
    st.download_button(
        label="Download JSON",
        data=json.dumps(json_data, indent=4),
        file_name='race_course_data.json',
        mime='application/json'
    )

    if st.button("Upload to Database"):
        if st.checkbox("Confirm ?"):
            coordinates = {
                "marks": [rc, pin, wg1, wg2],
                "boundary": boundary_coords[:-1]
            }
            insert_into_database(coordinates)
    

if __name__ == "__main__":
    main()
