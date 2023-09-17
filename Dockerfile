FROM public.ecr.aws/lambda/python:3.10

COPY lambda_function.py ${LAMBDA_TASK_ROOT}
COPY requirements.txt ${LAMBDA_TASK_ROOT}

RUN pip install --no-cache-dir -r requirements.txt

CMD ["lambda_function.lambda_handler"]