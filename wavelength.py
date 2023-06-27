import numpy as np
from matplotlib import pyplot as plt
import os
import math
import sys
import datetime
# Список имен файлов
fig = plt.figure("image")

sdate = sys.argv[1] #%d-%m-%Y
time = sys.argv[2] #%H:%M
filterid = sys.argv[3]
sdated = datetime.datetime.strptime(sdate + ' ' + time,'%d-%m-%Y %H:%M')
sdatet = sdated.strftime('%y%m%d')

# Получаем список имен файлов в заданной папке
path = "C:/Users/agey1/OneDrive/Рабочий стол" +"/" + sdatet
files = []
fdates = []

for f in os.listdir(path):
    if f.startswith("CE"+filterid+"_"+sdatet)  and f.endswith(".TIF"):
        files.append(f)

        fdates.append(datetime.datetime.strptime(f[4:16],'%y%m%d%H%M%S'))


images = []
for fnum in range(len(files)):
    if fdates[fnum] >= sdated:
        file = files[fnum]
        image = plt.imread(os.path.join(path, file))
        images.append(image)
print(images)
sum_images = np.zeros(images[0].shape)

for image in images:
    sum_images += image

group1 = images[:2]
group2 = images[3:5]
group3 = images[6:8]

# Считаем суммы изображений
sum1 = np.sum(group1, axis=0) / len(group1)
sum2 = np.sum(group2, axis=0) / len(group2)
sum3 = np.sum(group3, axis=0) / len(group3)

# Считаем разности изображений
diff1 = sum2 - sum1
diff2 = sum3 - sum2
diff = diff2 - diff1
#ce3 = r"C:\Users\agey1\OneDrive\Рабочий стол\221215\CE3_221215120120_0001.TIF"

# Вычисляем 2D- преобразование Фурье суммы изображений
fft_sum_images = np.fft.fft2(sum3)

# Вычисляем логарифмированный амплитудный спектр
spectrum = np.log(np.abs(np.fft.fftshift(fft_sum_images)))

# Рисуем изображения
figline = plt.figure("Intencity through line")
axline = figline.add_subplot(1, 1, 1)
fig, axs = plt.subplots(2, 3, figsize=(10, 8))
axs[0, 0].pcolormesh(sum1,vmin=0, vmax=30)
axs[0, 0].set_title("sum1")
axs[0, 1].pcolormesh(sum2)
axs[0, 1].set_title("sum2")
axs[0, 2].pcolormesh(sum3)
axs[0, 2].set_title("sum3")
axs[1, 0].pcolormesh(diff1, vmin=0, vmax=50)
axs[1, 0].set_title("diff1")
axs[1, 1].pcolormesh(diff2, vmin=0, vmax=50)
axs[1, 1].set_title("diff2")

## Функция для расчета экспоненциальной бегущей средней
def exp_running_avg(data, alpha):
  r_avg = [data[0]]
  for i in range(1, len(data)):
    r_avg.append(alpha * data[i] + (1-alpha) * r_avg[i-1])
  return np.array(r_avg)


# Выводим логарифмированный амплитудный спектр
axs[1, 2].pcolormesh(diff, cmap="gray", vmin=0, vmax=50)
axs[1, 2].set_title("diff")

# Объявляем переменные для хранения текущих координат точек и списка уже нарисованных отрезков
current_point = None
previous_point = None
lines = []

def onclick(event):
    print(event)
    global current_point, previous_point, lines
    if event.button == 1 and event.xdata is not None and event.ydata is not None:
        if current_point is not None:
            previous_point = current_point
        current_point = (int(round(event.xdata)), int(round(event.ydata)))
        if previous_point is not None:
            line = plt.plot([previous_point[0], current_point[0]], [previous_point[1], current_point[1]], 'r-')
            lines.append(line)
            fig.canvas.draw_idle()
            print(current_point, previous_point)

            y1 = previous_point[1]
            x1 = previous_point[0]

            y2 = current_point[1]
            x2 = current_point[0]

            a = (y2 - y1)/(x2-x1)
            b = y1 - a*x1
            if x1 < x2:
                x = np.arange(x1,x2,0.1)
            else:
                x = np.arange(x2,x1,0.1)
            print(x)
            y = (x*a + b)
            print(y)
            I = []
            l = []
            for i in range(len(x)):
                l.append(math.pow((x[i]-x[0])*(x[i]-x[0]) + (y[i]-y[0])*(y[i]-y[0]),0.5))
                I.append(diff[round(y[i]),round(x[i])])
            # Применяем к I экспоненциальную бегущую среднюю с коэффициентом alpha = 0.1
            I = exp_running_avg(I, 0.1)

            h = 100

            l2 = math.sqrt(math.pow((x2-x1),2) + math.pow((y2-y1),2))
            fi = l2*180/256

            L = 2*h*math.tan(fi*math.pi/180)

            print(L)

            print(l,I)
            axline.clear()
            distance = np.linspace(0,L,len(I))
            axline.plot(distance, I, marker="*", color="red")
            figline.canvas.draw()



    elif event.button == 3:
        if len(lines) > 0:
            previous_line = lines.pop()
            previous_line[0].remove()
            fig.canvas.draw_idle()

    return current_point, previous_point

# Связываем функцию-обработчик с событием нажатия кнопки мыши
fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()
