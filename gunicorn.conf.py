from gunicorn.glogging import Logger

bind = "0.0.0.0:5000"
wsgi_app = "app:create_app()"
accesslog = "-"
access_log_format = "%(h)s - - [%(t)s] %(r)s %(s)s %(b)s"


class HealthCheckFilteringLogger(Logger):
    # The Kubernetes readiness probe calls /health every few seconds, which buries the
    # requests worth reading. Successful probes are dropped; a failing one is still logged,
    # since that is the line you want when a pod drops out of the service.
    def access(self, resp, req, environ, request_time):
        if req.path == "/health" and resp.status_code == 200:
            return
        super().access(resp, req, environ, request_time)


logger_class = HealthCheckFilteringLogger
