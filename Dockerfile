FROM public.ecr.aws/python:3.10

COPY lambda_function.py {$LAMBDA_TASK_ROOT}
COPY requirements.txt {$LAMBDA_TASK_ROOT}

CMD ["python", "lambda_function.py"]