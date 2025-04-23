#!/home/eberendeev/anaconda3/bin/python3
#!/usr/bin/python3
#!/home/fano.icmmg/berendeev_e_a/anaconda3/bin/python3
#!/mnt/storage/home/eaberendeev/anaconda3/bin/python3
#!/mnt/storage/home/aaefimova/anaconda3/bin/python3

#from mpi4py import MPI

import time
import multiprocessing
import numpy as np
import os
import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as col
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
from matplotlib import rc
import pprint
import collections
import re

# Чтение параметров расчёта
from ReadDataLib import *
from LibPlot import *

# Конфигурационные параметры
CONFIG = {
    'work_dir': '../',
    'fields_amp': 0.03,  # амплитуда палитры для 3D полей
    'fields_amp_2d': 0.1,  # амплитуда палитры для 2D полей
    'dens_amp': 1.0,
    'bz0': 0.0,
    'step_y': "003",
    'i_start': 0,
    'i_max': 2999,
    'i_step': 50,
    'output_dir': '../Anime/2D/'
}

def create_output_directory(directory):
    """Создает директорию для вывода графиков, если она не существует"""
    try:
        os.makedirs(directory, exist_ok=True)
        print(f"Директория {directory} создана или уже существует")
    except OSError as e:
        print(f"Ошибка при создании директории {directory}: {e}")

def read_field_data(field_titles, full_base_name, time_step):
    """Чтение данных полей из файлов"""
    try:
        return ReadFieldsFile2DNew(full_base_name, field_titles, time_step)
    except Exception as e:
        print(f"Ошибка при чтении данных поля {full_base_name}: {e}")
        # Возвращаем пустой словарь вместо None
        exit(0)
        return {title: np.zeros((10, 10)) for title in field_titles}

def read_density_data(work_dir, particle_types, step_y, time_step):
    """Чтение данных плотности частиц"""
    density_data = collections.OrderedDict()
    for sort in particle_types:
        try:
            name = work_dir + "Particles/" + sort + "/Diag2D/Density_PlaneY_" + step_y + "_"
            density_data[sort] = ReadFieldsFile2DNew(name, ["Dens"], time_step)
        except Exception as e:
            print(f"Ошибка при чтении данных плотности для {sort}: {e}")
            density_data[sort] = np.zeros((10, 10))  # Default empty array
            exit(0)
    return density_data

def create_1d_slice(data_2d, axis='x', position=None):
    """Создает 1D срез из 2D данных"""
    if position is None:
        # Если позиция не указана, берем срез посередине
        if axis == 'x':
            position = data_2d.shape[1] // 2
            return data_2d[:, position]
        else:  # axis == 'z'
            position = data_2d.shape[0] // 2
            return data_2d[position, :]
    else:
        if axis == 'x':
            return data_2d[:, position]
        else:  # axis == 'z'
            return data_2d[position, :]

def plot_1d_slice(ax, data_1d, axis_label, title, color='blue'):
    """Строит 1D график среза"""
    x = np.arange(len(data_1d))
    ax.plot(x, data_1d, color=color)
    ax.set_xlabel(axis_label)
    ax.set_title(title)
    ax.grid(True)

def write_density_data(output_path, system_parameters, dens):
    """
    Записывает данные о плотности частиц в файл.
    
    Args:
        output_path (str): Путь к файлу для записи данных
        system_parameters (dict): Параметры системы с информацией о частицах
        dens (dict): Словарь со значениями плотности для каждого типа частиц
    """
    # Определяем, первый ли это запуск, проверяя существование файла
    first_run = not os.path.exists(output_path)
    
    try:
        # Создаем директорию, если она не существует
      #  os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        mode = "w" if first_run else "a"
        
        with open(output_path, mode) as file:
            if first_run:
                # Записываем заголовок с типами частиц при первом запуске
                file.write("# Particle density data\n")
                file.write("# " + "\t".join(system_parameters['Particles'].keys()) + "\n")
            
            # Записываем значения плотности (разделенные табуляцией)
            density_values = []
            for sort in system_parameters['Particles'].keys():
                density_values.append(str(dens[sort]))
            file.write("\t".join(density_values) + "\n")
            
        print(f"Density data written to {output_path}")
    except Exception as e:
        print(f"Error writing density data: {e}")

def main():
    # Загрузка системных параметров
    system_parameters = ReadParameters()
    system_parameters['2D_Diag_Fields'] = ['1', '1', '1', '1', '1', '1']
    
    # Создание директории для вывода
    create_output_directory(CONFIG['output_dir'])
    
    # Расчет временной задержки
    tdelay = int(system_parameters['TimeDelay2D']) * float(system_parameters['Dt'])
    
    print("Системные параметры:")
    print(system_parameters)
    
    # Основной цикл по временным шагам
    for i in range(CONFIG['i_start'], CONFIG['i_max'] + 1, CONFIG['i_step']):
        time_step = str(i).zfill(4)
        graph_name = 'ResXZ' + time_step + '.png'
        
        print('TimeStep=', time_step)
        
        # Создание фигуры с увеличенным размером для дополнительных графиков
        fig = plt.figure(figsize=(int(system_parameters['NumOfPartSpecies'])*4+40, 42))
        
        # Создание сетки для графиков с дополнительными столбцами для 1D графиков
        gs = gridspec.GridSpec(6, 4, height_ratios=[0.1, 1, 1, 1, 1, 1])
        
        # Функция для создания подграфика
        def sub_plot(x, y):
            return fig.add_subplot(gs[y, x])
        
        # Чтение данных полей
        field_e_base = CONFIG['work_dir'] + "Fields/Diag2D/FieldE_PlaneY_" + CONFIG['step_y'] + "_"
        field_e = read_field_data(["Ex", "Ey", "Ez"], field_e_base, time_step)
        
        field_b_base = CONFIG['work_dir'] + "Fields/Diag2D/FieldB_PlaneY_" + CONFIG['step_y'] + "_"
        field_b = read_field_data(["Bx", "By", "Bz"], field_b_base, time_step)
        
        # Чтение данных плотности
        dens_xy = read_density_data(CONFIG['work_dir'], system_parameters['Particles'].keys(), CONFIG['step_y'], time_step)
        
        # Построение 2D графиков полей
        Plot2Ddens(field_e["Ex"], "xz", -CONFIG['fields_amp'], CONFIG['fields_amp'], 
                  "Ex", "field", 1, system_parameters, sub_plot(2, 1), fig)
        # Plot2Ddens(field_e["Ey"], "xz", -CONFIG['fields_amp'], CONFIG['fields_amp'], 
        #           "Ey", "field", 1, system_parameters, sub_plot(2, 1), fig)
        Plot2Ddens(field_e["Ez"], "xz", -CONFIG['fields_amp'], CONFIG['fields_amp'], 
                  "Ey", "field", 1, system_parameters, sub_plot(2, 2), fig)
        Plot2Ddens(field_b["Bx"], "xz", -CONFIG['fields_amp']+CONFIG['bz0'], CONFIG['fields_amp']+CONFIG['bz0'], 
                  "Bx", "field", 1, system_parameters, sub_plot(2, 3), fig)
        Plot2Ddens(field_b["Bz"], "xz", -CONFIG['fields_amp']+CONFIG['bz0'], CONFIG['fields_amp']+CONFIG['bz0'], 
                  "Bz", "field", 1, system_parameters, sub_plot(2, 4), fig)
        
        # Добавление 1D графиков для полей (срезы по центру)
        # Срезы по X
        #ex_slice_x = create_1d_slice(field_e["Ex"], axis='x')
        #ey_slice_x = create_1d_slice(field_e["Ey"], axis='x')
       # bx_slice_x = create_1d_slice(field_b["Bx"], axis='x')
       # bz_slice_x = create_1d_slice(field_b["Bz"], axis='x')
        
        #plot_1d_slice(sub_plot(1, 1), ex_slice_x, "Z", "Ex slice along Z", color='red')
        #plot_1d_slice(sub_plot(3, 1), ey_slice_x, "Z", "Ey slice along Z", color='red')
        #plot_1d_slice(sub_plot(5, 1), bx_slice_x, "Z", "Bx slice along Z", color='red')
        #plot_1d_slice(sub_plot(7, 1), bz_slice_x, "Z", "Bz slice along Z", color='red')
        
        # Срезы по Z
        ex_slice_z = create_1d_slice(field_e["Ex"], axis='z')
        ey_slice_z = create_1d_slice(field_e["Ey"], axis='z')
        ez_slice_z = create_1d_slice(field_e["Ey"], axis='z')
        bx_slice_z = create_1d_slice(field_b["Bx"], axis='z')
        bz_slice_z = create_1d_slice(field_b["Bz"], axis='z')
        
        plot_1d_slice(sub_plot(3, 1), ex_slice_z, "X", "Ex slice along X", color='blue')
        #plot_1d_slice(sub_plot(3, 1), ey_slice_z, "X", "Ey slice along X", color='blue')
        plot_1d_slice(sub_plot(3, 2), ez_slice_z, "X", "Ez slice along X", color='blue')
        plot_1d_slice(sub_plot(3, 3), bx_slice_z, "X", "Bx slice along X", color='blue')
        plot_1d_slice(sub_plot(3, 4), bz_slice_z, "X", "Bz slice along X", color='blue')
        
        # Построение графиков плотности частиц
        gs_y = 1
        dens = {}
        for sort in system_parameters['Particles'].keys():
            # 2D график плотности
            Plot2Ddens(np.abs(dens_xy[sort]), "xz", 0, 
                      CONFIG['dens_amp'] * float(system_parameters['Particles'][sort]['Dens']),
                      sort + " density", "dens", 1, system_parameters, sub_plot(0, gs_y), fig)
            sort_dens = np.abs(dens_xy[sort])[2:-3]
            dens[sort] = np.mean(sort_dens)
            # 1D срезы плотности
            dens_slice_x = create_1d_slice(np.abs(dens_xy[sort]), axis='x')
            dens_slice_z = create_1d_slice(np.abs(dens_xy[sort]), axis='z')
            
            #plot_1d_slice(sub_plot(1, gs_y), dens_slice_x, "Z", f"{sort} density slice along Z", color='red')
            plot_1d_slice(sub_plot(1, gs_y), dens_slice_z, "X", f"{sort} density slice along X", color='blue')
            
            gs_y += 1
        
        # Добавление информации о времени
        max_time = int(time_step) * tdelay
        max_time_dt = round(max_time, 1)
        
        plt.figtext(0.5, 0.96, r'$t\cdot\omega_p = ' + str(max_time_dt) + '$', 
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                   color='black', fontsize=18, ha='center')
        
        plt.tight_layout()
        plt.savefig(CONFIG['output_dir'] + graph_name, format='png', dpi=150)
        plt.close(fig)
        # Запись данных о плотности в файл
        write_density_data(CONFIG['output_dir'] + "dens_data.txt", system_parameters, dens)
        

if __name__ == "__main__":
    main()