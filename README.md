# Репозиторій з Лабораторними роботами з премету "Теоретико-числові алгоритми в криптології"

**Варіант:** 2

**Виконали:** Бондар Петро, Кістаєв Матвій

**Група:** ФІ-03

## Lab 1: Пошук канонічного розкладу великого числа, використовуючи відомі методи факторизації

1. Для інтерактивної і більш зручної взаємодії радимо використовувати jupyter-notebook.
Відкрийте файл `number-factorization.ipynb` у будь-якому зручному редакторі, який підтримує такий формат.

2. Для запуску програми за допомогою консолі:

        python script.py [<prime base limit>] [<sieving range (M)>] [<max smooth numbers>]

3. Запуск за допомогою Docker (https://hub.docker.com/r/petrob2003/nta-bondar-kistaiev-lab-1):
        
        # Pull image
        docker pull petrob2003/nta-bondar-kistaiev-lab-1

        # run image (add parameters if needed)
        docker run -i petrob2003/nta-bondar-kistaiev-lab-1 [<prime base limit>] [<sieving range (M)>] [<max smooth numbers>]