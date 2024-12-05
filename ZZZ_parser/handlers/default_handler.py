from ZZZ_parser.handlers.base_handler import BaseHandler


class DefaultHandler(BaseHandler):
    def can_handle(self, data):
        print("placeholder default handler")

    def handle(self, data):
        print("placeholder default handler")
