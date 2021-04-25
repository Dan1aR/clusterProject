# clusterProject
Проект по кластеризации, Андреев Данил ИУ6-23Б, МГТУ им. Н.Э.Баумана. 
Алгоритмы: K-средних, Иерархический, Волновой

Описание:

main.cpp - программа, осуществляющая кластеризацию
-
Содержит интерфейс для:
Добавления точек в виде "облаков" (Эллипсов)
Добавления точек из файла
Запуск Алгоритмов поиска кластеров

Результат сохраняется в output.txt по такому же шаблону

analytics.ipynb - Удобная визуализация

# gnuplot
Рекомендую следующую команду для визуализации:
plot "< awk '{if($3 == \"0\") print}' output.txt" u 1:2 t "red" w p pt 2, "< awk '{if($3 == \"1\") print}' output.txt" u 1:2 t "green" w p pt 2, "< awk '{if($3 == \"2\") print}' output.txt" u 1:2 t "blue" w p pt 2, "< awk '{if($3 == \"3\") print}' output.txt" u 1:2 t "pink" w p pt 2
